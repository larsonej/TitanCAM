      subroutine lognormal( totn ,rn, rsig )
c
c  
c  @(#) lognorm.f  Ackerman  Sep-1997
c
c  This routine fits a log-normal distribution to the particle number
c  distributions for each particle group at each grid point.
c
c
c  Argument list input:
c    None.
c
c  Argument list output:
c
c    totn    Total particle concentration
c    rn      Number mode radius (geometric mean radius by number)
c    rsig    Geometric standard deviation
c
c
c  Include global constants and variables
c  
      include 'globaer.h'
c
c
c  Declare local variables
c
      dimension totn(NXYZ,NGROUP), rn(NXYZ,NGROUP), rsig(NXYZ,NGROUP),
     $          pc_cgs(NBIN)
c
c
c-------------------------------------------------------------------------------
c  
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter lognorm'
c
c-------------------------------------------------------------------------------
c
c
c
c  Fit log-normal distribution to particle concentration distributions
c  at each grid point.
c
      do ig = 1,NGROUP
        ie = ienconc(ig)
        do ixyz = 1,NXYZ
c
c
c  Convert concentrations to cgs units and compute total
c  particle concentration, <totn>, and total of n*log(r), <totlnr>
c
          totn(ixyz,ig) = 0.
          totlnr = 0.

          do i = 1,NBIN

            pc_cgs(i) = pc3(ixyz,i,ie) /
     $        ( xmet3(ixyz)*ymet3(ixyz)*zmet3(ixyz) )

            totn(ixyz,ig) = totn(ixyz,ig) + pc_cgs(i)

            totlnr = totlnr + pc_cgs(i) * log( r(i,ig) )

          enddo
c
c
c  Compute geometric mean number radius <rn>
c  and geometric standard deviation <rsig>
c
          if( totn(ixyz,ig) .gt. (NBIN+1)*SMALL_PC )then

             alnrn = totlnr / totn(ixyz,ig)
             rn(ixyz,ig) = exp( alnrn )

             totlnrsig = 0.
             do i = 1,NBIN
              totlnrsig = totlnrsig + 
     $                    pc_cgs(i) * ( log( r(i,ig) ) - alnrn )**2
             enddo

             rsig(ixyz,ig) = exp( sqrt( totlnrsig / totn(ixyz,ig) ) )

          else

             rn(ixyz,ig) = 0.
             rsig(ixyz,ig) = 0.

          endif

        enddo    

      enddo    

      return
      end
