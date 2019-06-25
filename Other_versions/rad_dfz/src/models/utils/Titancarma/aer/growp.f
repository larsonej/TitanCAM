      subroutine growp(ibin,ielem)
c
c
c  @(#) growp.f  Ackerman  Dec-1995
c  This routine calculates particle source terms due to growth <growpe>
c  for one particle size bin at one spatial grid point per call.
c  Globals <ix>, <iy>, <iz>, <ixy>, <ixyz> specify current spatial pt's indices.
c
c
c  Argument list input:
c    ibin, ielem
c
c  Argument list output:
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Local declarations
c
c
c  Define formats
c
    2 format(a,i4,3x,'growlg: ',1pe11.4,3x,'add: ',1pe11.4,3x,
     $       'dmdt_gc: ',1pe11.4,3x,'dmdt_pr: ',1pe11.4,i10)
    3 format(a,i4,3x,'gprod_grow: ',1pe11.4,'  add: ',1pe11.4,
     $       3x,'bin:',i3)
    4 format(a,'  from vapor(add): ',1pe11.4,'  total(growcore): ',
     $       1pe11.4,'  total(primary): ',1pe11.4,
     $       '  gprod_grow(gc): ',1pe11.4)
    5 format(a,i4,3x,'growlg: ',1pe11.4,3x,'add: ',1pe11.4,i10)
    6 format(a,i4,3x,'both: ',1pe11.4,3x,'prim: ',1pe11.4,3x,
     $       'gcor:',1pe11.4)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter growp'
c

c  Define group & particle # concentration indices for current element
c
      igroup = igelem(ielem)      ! target particle group 
      iepart = ienconc(igroup)	  ! target particle number concentration element
c
c
c  Calculate production terms due to condensational growth <growpe>
c  only if group to which element belongs grows.
c
      if( igrowgas(iepart) .ne. 0 .and. ibin .ne. 1 )then
c
c
c  Bypass calculation if few droplets are present 
c
         if( pconmax(ixyz,igroup) .gt. FEW_PC )then

          growpe(ibin,ielem) = pc3(ixyz,ibin-1,ielem)
     $                            * growlg(ibin-1,igroup) 
          
         if( itype(ielem) .eq. I_VOLCORE ) then  
           igas = igrowgas(ielem-1) 
           if( t3(ixyz) .le. Tfreez(igas)) 
     $      growpe(ibin,ielem) = growpe(ibin,ielem) +
     $           pc3(ixyz,ibin-1,iepart)*growlg(ibin-1,igroup)* 
     $                  diffmass(ibin,igroup,ibin-1,igroup) 
         endif
          
          if( itype(ielem) .eq. I_GROWCORE ) then

           igas = igrowgas(ielem)


c     <add> is g/cm3 of growcore gas which is removed from vapor phase
c     (and added to particle in ibin) due to growth so that the fraction
c     will stay constant as the particle moves across bins

           dmdt_main = dmdt_gro(ibin-1,igroup) 
     $                     - dmdte_gro(ibin-1,ielem)
           if( dmdt_main  .gt. ZERO ) then
            if( dmdte_gro(ibin-1,ielem) .gt. ZERO ) then

             ! Maingas and growcore gas should both grow on particle
             ! (for now will grow at constant fraction) 

             add = pc3(ixyz,ibin-1,iepart) * growlg(ibin-1,igroup) *
     $               diffmass(ibin,igroup,ibin-1,igroup) *
     $               dmdte_gro(ibin-1,ielem)/dmdt_gro(ibin-1,igroup)


c            if(ixyz.eq.iwa)
c    $        write(*,5) 'Both grow',
c    $                   ibin-1,growlg(ibin-1,igroup),add,itime

            else
  
             ! Maingas should grow, but growcore gas should evaporate

c            if(ixyz.eq.iwa)
c    $         write(*,*) 'Grow only maingas',
c    $                     ibin-1,growlg(ibin-1,igroup),itime
             add = ZERO

            endif

           else
 
            ! Grow the growcore, but keep the primary cloud volatile 
            ! constant
            ! (might need to change this part for multiple growcores, ie.
            !  add condition if current growcore should not grow)

            add = pc3(ixyz,ibin-1,iepart)*growlg(ibin-1,igroup)* 
     $                  diffmass(ibin,igroup,ibin-1,igroup) 
 
c            if(ixyz.eq.iwa)
c    $         write(*,2) 'Growcore only grows',ibin-1,
c    $          growlg(ibin-1,igroup),add,dmdte_gro(ibin-1,ielem),
c    $          dmdt_main,itime

           endif

           growpe(ibin,ielem) = growpe(ibin,ielem) + add
           gprod_grow(igroup,igas) = gprod_grow(igroup,igas) - add

c          if(ixyz.eq.iwa)
c    $       write(*,3) '<growp> gprod_grow',igas,
c    $                  gprod_grow(igroup,igas),add,ibin

          endif !GROWCORE
         endif !.gt. FEW_PC
      endif
c
c
c  Return to caller with growth production terms evaluated.
c
      return
      end
