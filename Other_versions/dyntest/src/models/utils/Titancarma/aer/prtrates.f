      subroutine prtrates
c
c
c  @(#) prtrates.f  Barth  Jan-2001
c  This routine outputs particle production and loss rates at the current 
c  timestep to an output print file.
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define formats
c
c   1 format(f10.2,f7.2,2(i3),i4,1pe10.2)
    1 format(f12.6,f5.1,2(i3),i4,1pe10.2)
c
c
       character*(50) ratefile
c
c
c---------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter prtrates'
      return
c
c
      if(ixyz .eq. 1) then
c
c   ..Initial file names:
c
c   ..First file name tag to be appended to initial file name
        iun=iun+1
        if (iun.gt.9) then
          iun=0
c   ..Second file name tag
          ideux=ideux+1
          if (ideux.gt.9) then
            ideux=0
            itrois=itrois+1
            if (itrois.gt.9) then
              itrois=0
              iquatre=iquatre+1
            endif
          endif
        endif
c
c  ..Construct file names..
c
        ratefile='Files/Output/rates'//ext//char(iquatre+48)//
     $   char(itrois+48)//char(ideux+48)//char(iun+48)//'.p'
c
c   ..Report file name writen to...
c
        print*,"Writing to:"
        print*,ratefile
c
c
c  ..Write out production and loss rates
c
        open(unit=2,file=ratefile,STATUS='UNKNOWN')
       
      endif

       do ielem = 1,2  !only write out rates for tholin and ice
           do ibin = 1,NBIN

              igroup = igelem(ielem)         ! particle group

              rnuclgtot = 0.
              do igto = 1,NGROUP
                  rnuclgtot = rnuclgtot + rnuclg(ibin,igroup,igto)
              enddo

c       rates in cm-3 s-1 
c
c
c        Coagulation (production only, id = 1)
c
            if( coagpe(ixyz,ibin,ielem) .ge. SMALL_PC )
     $        write(2,1) time/60.**2/24.,zl3(ixyz)/1.d5,1,
     $           ielem,ibin,coagpe(ixyz,ibin,ielem)
c
c        Nucleation (tholin loss only, id = 2)
c
            if( ielem .eq. 1 .and. rnuclgtot .ne. 0. ) 
     $        write(2,1) time/60.**2/24.,zl3(ixyz)/1.d5,2,
     $          ielem,ibin,rnuclgtot*pc3(ixyz,ibin,ielem)
c
c        Growth (cloud production only, id = 3)
c
            if( ielem .eq. 2 .and. growpe(ibin,ielem) .ge. SMALL_PC )
     $      write(2,1) time/60.**2/24.,zl3(ixyz)/1.d5,3,
     $        ielem,ibin,growpe(ibin,ielem)
c             
c        Evaporation (tholin production id = 4, cloud loss id = 5)
c
              if(ielem .eq. 1 .and. evappe(ibin,ielem) .ge. SMALL_PC)
     $        write(2,1) time/60.**2/24.,zl3(ixyz)/1.d5,4,
     $          ielem,ibin,evappe(ibin,ielem)
              if(ielem .eq. 2 .and. evaplg(ibin,igroup) .ge. SMALL_PC)
     $        write(2,1) time/60.**2/24.,zl3(ixyz)/1.d5,5,
     $          ielem,ibin,evaplg(ibin,igroup)*pc3(ixyz,ibin,ielem)

          enddo
      enddo

      if(ixyz .eq. NXYZ) close(unit=2)
c
c
c  Return to caller with rates info output to print file
c
      return
      end
