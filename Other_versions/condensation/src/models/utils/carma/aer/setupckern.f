       subroutine setupckern 
c 
c 
c  @(#) setupckern.f  Ackerman Oct-1995 
c 
c  This routine evaluates the coagulation kernels, ckernel(k,j1,j2,i1,i2)
c  [cm^-3 s^-1]. Indices correspond to vertical level <k>, aerosol groups
c  <j1,j2> and bins <i1,i2> of colliding particles.
c
c  This routine requires that vertical profiles of temperature <T>,
c  air density <rhoa>, and viscosity <rmu> are defined.
c  (i.e., initatm.f must be called before this)
c  The vertical profile with ix = iy = 1 is used.
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
c  Local declarations
c
c    <e_coll2> is 2-D collision efficiency for current group pair under
c    consideration (for extrapolation of input data)
c
      dimension e_coll2(NBIN,NBIN)
c
c     
c    <fill_bot> is used for filling <icoag>
c
      logical fill_bot
c
c     
c    <e_small> is smallest collision efficiency we use when interpolating data

      parameter ( e_small = 0.0001 )
c     
c
c    <NP_DATA> is number of collector/collected pairs in input data 
c    <NR_DATA> is number of radius bins in input data

      parameter ( NP_DATA = 21, NR_DATA = 12 )  !EJL - should be changed to match number of radii bins?
c
c    <data_p> are radius ratios (collected/collector)
c    <data_r> are collector drop radii (um)
c    <data_e> are geometric collection efficiencies
c
      dimension data_p(NP_DATA), data_r(NR_DATA), 
     $          data_e(NP_DATA,NR_DATA)
c
c
c  Initialization of input data for gravitational collection.
c  The data were compiled by Hall (J. Atmos. Sci. 37, 2486-2507, 1980).
c
      data data_p/0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,
     1 0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00/
c
      data data_r( 1), (data_e(ip, 1),ip=1,NP_DATA) /   10.0,
     $ 0.0001, 0.0001, 0.0001, 0.0001, 0.0140, 0.0170, 0.0190, 0.0220, 
     $ 0.0270, 0.0300, 0.0330, 0.0350, 0.0370, 0.0380, 0.0380, 0.0370, 
     $ 0.0360, 0.0350, 0.0320, 0.0290, 0.0270 /
      data data_r( 2), (data_e(ip, 2),ip=1,NP_DATA) /   20.0,
     $ 0.0001, 0.0001, 0.0001, 0.0050, 0.0160, 0.0220, 0.0300, 0.0430, 
     $ 0.0520, 0.0640, 0.0720, 0.0790, 0.0820, 0.0800, 0.0760, 0.0670, 
     $ 0.0570, 0.0480, 0.0400, 0.0330, 0.0270 /
      data data_r( 3), (data_e(ip, 3),ip=1,NP_DATA) /   30.0,
     $ 0.0001, 0.0001, 0.0020, 0.0200, 0.0400, 0.0850, 0.1700, 0.2700, 
     $ 0.4000, 0.5000, 0.5500, 0.5800, 0.5900, 0.5800, 0.5400, 0.5100, 
     $ 0.4900, 0.4700, 0.4500, 0.4700, 0.5200 /
      data data_r( 4), (data_e(ip, 4),ip=1,NP_DATA) /   40.0,
     $ 0.0001, 0.0010, 0.0700, 0.2800, 0.5000, 0.6200, 0.6800, 0.7400, 
     $ 0.7800, 0.8000, 0.8000, 0.8000, 0.7800, 0.7700, 0.7600, 0.7700, 
     $ 0.7700, 0.7800, 0.7900, 0.9500, 1.4000 /
      data data_r( 5), (data_e(ip, 5),ip=1,NP_DATA) /   50.0,
     $ 0.0001, 0.0050, 0.4000, 0.6000, 0.7000, 0.7800, 0.8300, 0.8600, 
     $ 0.8800, 0.9000, 0.9000, 0.9000, 0.9000, 0.8900, 0.8800, 0.8800, 
     $ 0.8900, 0.9200, 1.0100, 1.3000, 2.3000 /
      data data_r( 6), (data_e(ip, 6),ip=1,NP_DATA) /   60.0,
     $ 0.0001, 0.0500, 0.4300, 0.6400, 0.7700, 0.8400, 0.8700, 0.8900, 
     $ 0.9000, 0.9100, 0.9100, 0.9100, 0.9100, 0.9100, 0.9200, 0.9300, 
     $ 0.9500, 1.0000, 1.0300, 1.7000, 3.0000 /
      data data_r( 7), (data_e(ip, 7),ip=1,NP_DATA) /   70.0,
     $ 0.0001, 0.2000, 0.5800, 0.7500, 0.8400, 0.8800, 0.9000, 0.9200, 
     $ 0.9400, 0.9500, 0.9500, 0.9500, 0.9500, 0.9500, 0.9500, 0.9700, 
     $ 1.0000, 1.0200, 1.0400, 2.3000, 4.0000 /
      data data_r( 8), (data_e(ip, 8),ip=1,NP_DATA) /  100.0,
     $ 0.0001, 0.5000, 0.7900, 0.9100, 0.9500, 0.9500, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      data data_r( 9), (data_e(ip, 9),ip=1,NP_DATA) /  150.0,
     $ 0.0001, 0.7700, 0.9300, 0.9700, 0.9700, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      data data_r(10), (data_e(ip,10),ip=1,NP_DATA) /  200.0,
     $ 0.0001, 0.8700, 0.9600, 0.9800, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      data data_r(11), (data_e(ip,11),ip=1,NP_DATA) /  300.0,
     $ 0.0001, 0.9700, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
      data data_r(12), (data_e(ip,12),ip=1,NP_DATA) / 1000.0,
     $ 0.0001, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 
     $ 1.0000, 1.0000, 1.0000, 1.0000, 1.0000 /
c
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupckern'
c
c-------------------------------------------------------------------------------

c
c      electron charge [esu]
        e = -4.8d-10
c
c
c  Fill <icoag>, maintaining diagonal symmetry
c
c  Fill bottom of matrix if non-zero term(s) in upper half;
c  also check for non-zero, non-matching, non-diagonal terms.
c
      fill_bot = .true.
      do irow = 2, NGROUP
        do icol = 1, irow-1
          if( icoag(irow,icol) .ne. 0 )then
            fill_bot = .false.
            if( icoag(icol,irow) .ne. 0 .and.
     $          icoag(icol,irow) .ne. icoag(irow,icol) )then
              print*,'stop in setupckern(): bad icoag array'
              call endcarma
            endif
          endif
        enddo
      enddo

      do ig = 2, NGROUP
        do jg = 1, ig-1
          if( fill_bot )then
            irow = ig
            icol = jg
          else
            irow = jg
            icol = ig
          endif
          icoag(irow,icol) = icoag(icol,irow)
        enddo
      enddo
c
c
c  <cstick> is the probability that two particles that collide
c  through thermal coagulation will stick to each other.
c
!      cstick = 1.  !EJL
c
c
c  Use constant kernel if <icoagop> = 0
c
      if( icoagop .eq. 0 )then
        do j1 = 1, NGROUP
         do j2 = 1, NGROUP
          do i1 = 1, NBIN
           do i2 = 1, NBIN
            do k = 1, NZ
             ckernel(k,i1,i2,j1,j2) = ck0
            enddo
           enddo
          enddo
         enddo
        enddo
        return   ! Return to caller with coagulation kernels evaluated.
      endif
 
!      write(99,*) 'ckern1'
!	  write(99,*) ckernel  !EJL
 
      if( icollec .eq. 2 )then
c
c   Convert <data_r> from um to cm and take logarithm of <data_e>.
c
        do i = 1, NR_DATA
          data_r(i) = data_r(i)/1.e4
          do ip = 1, NP_DATA
            data_e(ip,i) = log(data_e(ip,i))
          enddo
        enddo

      endif
 
c
c  Loop over vertical layers (use column ix = iy = 1)
c
      ix = 1
      iy = 1

      do k = 1, NZ
c
c  This is <rhoa> in cartesian coordinates.
c
        rhoa_cgs = rhoa(ix,iy,k) /
     $             (xmet(ix,iy,k)*ymet(ix,iy,k)*zmet(ix,iy,k))

        temp1 = BK*t(ix,iy,k)
        temp2 = 6.*PI*rmu(k)

c  prrat is the charge to radius ratio. Lavvas et al. 2010 used 15.
        prrat(k) = 15*1.e4   !EJL 15 e/um * 1.e4 um/cm
c
c
c  Loop over groups
c
        do j1 = 1, NGROUP
        do j2 = 1, NGROUP

        if( icoag(j1,j2) .ne. 0 )then

c
c  First particle
c
         do i1 = 1, NBIN

          r1 = rf(i1,j1)
          di = temp1*bpm(k,i1,j1)/(temp2*r1)
          gi  = sqrt( 8.*temp1/(PI*rmass(i1,j1)) )
          rlbi = 8.*di/(PI*gi)
          dti1= (2.*r1 + rlbi)**3
          dti2= (4.*r1*r1 + rlbi*rlbi)**1.5
          dti = 1./(6.*r1*rlbi)
          dti = dti*(dti1 - dti2) - 2.*r1

          do i2 = 1, NBIN
c
c  Second particle
c
            r2  = rf(i2,j2)
            dj  = temp1*bpm(k,i2,j2)/(temp2*r2)
            gj  = sqrt( 8.*temp1/(PI*rmass(i2,j2)) )
            rlbj = 8.*dj/(PI*gj)
            dtj1= (2.*r2 + rlbj)**3
            dtj2= (4.*r2*r2 + rlbj*rlbj)**1.5
            dtj = 1./(6.*r2*rlbj)
            dtj = dtj*(dtj1 - dtj2) - 2.*r2
c
c  First calculate thermal coagulation kernel
c  
            rp  = r1 + r2
            dp  = di + dj
            gg0  = sqrt(gi*gi + gj*gj)
            delt= sqrt(dti*dti + dtj*dtj)
!            term1 = rp/(rp + delt)
!            term2 = 4.*dp/(gg*rp)
            cbr_term0(k,i1,i2,j1,j2) = 4.*PI*rp*dp
            cbr_term1(k,i1,i2,j1,j2) = rp/(rp + delt)
            cbr_term2(k,i1,i2,j1,j2) = 4.*dp/(gg0*rp)
c
c   EJL
c   <cstick> is probability that two particles that collide
c   through thermal coagulation will stick stick together.
             cstick = exp(-r1*r2*prrat(k)*prrat(k)*e*e/(temp1*rp))
c
c   <cbr> is thermal (brownian) coagulation coefficient
c
!            cbr = 4.*PI*rp*dp/(term1 + term2)

! From E.Barth
             cbr = cbr_term0(k,i1,i2,j1,j2) /
     $               (cbr_term1(k,i1,i2,j1,j2)
     $                  + cbr_term2(k,i1,i2,j1,j2)/cstick)
	 
	 !Write out cbr
!      open(unit=99, file='carma_ckern.txt', status='unknown') !EJL
!      write(99,*) cbr
!	  close(unit=99)	 
c 
c   Determine indices of larger and smaller particles (of the pair)
c
            if (r2 .ge. r1) then
              r_larg = r2
              r_smal = r1
              i_larg = i2
              i_smal = i1
              ig_larg = j2
              ig_smal = j1
            else
              r_larg = r1
              r_smal = r2
              i_larg = i1
              i_smal = i2
              ig_larg = j1
              ig_smal = j2
            endif
c 
c   Calculate enhancement of coagulation due to convective diffusion 
c   as described in Pruppacher and Klett.
c
c   Enhancement applies to larger particle.
c
            re_larg = re(k,i_larg,ig_larg)
c
c   <pe> is Peclet number.
c
            pe  = re_larg*rmu(k) / (rhoa_cgs*di)
            pe3 = pe**(1./3.)
c
c   <ccd> is convective diffusion coagulation coefficient
c
            if( re_larg .lt. 1. )then
              ccd = 0.45*cbr*pe3
            else 
              ccd = 0.45*cbr*pe3*re_larg**(ONE/6.)
            endif
c
c   Next calculate gravitational collection kernel.  
c
c   First evaluate collection efficiency <e>.
c 
            if( icollec .eq. 0 )then
c
c   constant value
c
              e_coll = grav_e_coll0

            else if( icollec .eq. 1 )then
c
c   Find maximum of Langmuir's formulation and Fuchs' value.
c   First calculate Langmuir's efficiency <e_langmuir>.
c
c   <sk> is stokes number.
c   <vfc_{larg,smal}> is the fallspeed in cartesian coordinates.
c
              vfc_smal = vf(k,i_smal,ig_smal) * zmet2(1,k)
              vfc_larg = vf(k,i_larg,ig_larg) * zmet2(1,k)

              sk = vfc_smal * (vfc_larg - vfc_smal) / (r_larg*GRAV)
 
              if( sk .lt. 0.08333334 )then
                e1 = 0.
              else 
                e1 = (sk/(sk + 0.25))**2
              endif
 
              if( sk .lt. 1.214 )then
                e3  = 0.
              else
                e3  = 1./(1.+.75*log(2.*sk)/(sk-1.214))**2
              endif
 
              if( re_larg .lt. 1. )then
                e_langmuir = e3
              else if( re_larg .gt. 1000. )then
                e_langmuir = e1
              else if( re_larg .le. 1000. )then
                re60 = re_larg/60.
                e_langmuir = (e3  + re60*e1)/(1. + re60)
              endif
c
c   Next calculate Fuchs' efficiency (valid for r < 10 um).
c
              pr = r_smal/r_larg
              e_fuchs   = (pr/(1.414*(1. + pr)))**2

              e_coll = max( e_fuchs, e_langmuir )
 
 
 
            else if( icollec .eq. 2 )then
c
c   Interpolate input data (from data statment at beginning of subroutine).
c
              pr = r_smal/r_larg
c
c
c   First treat cases outside the data range
c
              if( pr .lt. data_p(2) )then
c
c   Radius ratio is smaller than lowest nonzero ratio in input data --
c   use constant values (as in Beard and Ochs, 1984) if available,
c   otherwise use very small efficiencty
c
                if( i2 .eq. i_larg )then
                  if( i2.eq.1 )then
                    e_coll = e_small
                  else
                    e_coll = e_coll2(i1,i2-1)
                  endif
                else
c!  EJL changed i2 to i1 in line below				
                  if( i1.eq.1 )then  
                    e_coll = e_small
                  else
                    e_coll = e_coll2(i1-1,i2)
                  endif
                endif

              elseif( r_larg .lt. data_r(1) )then
c
c   Radius of larger particle is smaller than smallest radius in input data -- 
c   assign very small efficiency.
c
                  e_coll = e_small

                else
c
c
c   Both droplets are either within grid (interpolate) or larger
c   droplet is larger than maximum on grid (extrapolate) -- in both cases 
c   will interpolate on ratio of droplet radii.
c
c   Find <jp> such that data_p(jp) <= pr <= data_p(jp+1) and calculate
c   <pblni> = fractional distance of <pr> between points in <data_p> 
c
                jp = NP_DATA
                do jj = NP_DATA-1, 2, -1
                  if( pr .le. data_p(jj+1) ) jp = jj
                enddo
c should not need this if-stmt
                if( jp .lt. NP_DATA )then
                  pblni = (pr - data_p(jp)) 
     $                  / (data_p(jp+1) - data_p(jp))
                else
c nor this else-stmt
                  print*,'stop in setupckern: NP_DATA < jp = ',jp
                  call endcarma
                endif

                if( r_larg .gt. data_r(NR_DATA) )then
c
c    Extrapolate on R and interpolate on p 
c
                  e_coll = (1.-pblni)*data_e(jp  ,jr) +
     $                     (   pblni)*data_e(jp+1,jr)

                else
c
c    Find <jr> such that data_r(jr) <= r_larg <= data_r(jr+1) and calculate
c    <rblni> = fractional distance of <r_larg> between points in <data_r>
c
                  jr = NR_DATA
                  do jj = NR_DATA-1, 1, -1
                    if( r_larg .le. data_r(jj+1) ) jr = jj
                  enddo
                  rblni = (r_larg - data_r(jr)) 
     $                  / (data_r(jr+1) - data_r(jr))
c             
c    Bilinear interpolation of logarithm of data.
c
                  e_coll = (1.-pblni)*(1.-rblni)*data_e(jp  ,jr  ) +
     $                     (   pblni)*(1.-rblni)*data_e(jp+1,jr  ) +
     $                     (1.-pblni)*(   rblni)*data_e(jp  ,jr+1) +
     $                     (   pblni)*(   rblni)*data_e(jp+1,jr+1)

c    (since data_e is logarithm of efficiencies)

                  term1 = (1.-rblni)*(1.-pblni)*data_e(jp,jr)
                  if( jp .lt. NP_DATA )then
                    term2 = pblni*(1.-rblni)*data_e(jp+1,jr)
                  else
                    term2 = -100.
                  endif
                  if( jr .lt. NR_DATA )then
                    term3 = (1.-pblni)*rblni*data_e(jp,jr+1)
                  else
                    term3 = -100.
                  endif
                  if( jr .lt. NR_DATA .and. jp .lt. NP_DATA )then
                    term4 = pblni*rblni*data_e(jp+1,jr+1)
                  else
                    term4 = -100.
                  endif
    
                  e_coll = exp(term1 + term2 + term3 + term4)
    
                endif
              endif

              e_coll2(i1,i2) = e_coll

            endif
c
c  Now calculate coalescence efficiency from Beard and Ochs 
c  (J. Geophys. Res. 89, 7165-7169, 1984).
c  
            beta = log(r_smal*1.e4) + 0.44*log(r_larg*50.)
            b_coal = 0.0946*beta - 0.319
            a_coal = sqrt(b_coal**2 + 0.00441)
            x_coal = (a_coal-b_coal)**(ONE/3.) 
     $             - (a_coal+b_coal)**(ONE/3.)
            x_coal = x_coal + 0.459
c
c  Limit extrapolated values to no less than 50% and no more than 100%
c
            x_coal = max(x_coal,.5*ONE)
            e_coal = min(x_coal,1.*ONE)
c
c  Now use coalescence efficiency and collision efficiency in definition
c  of (geometric) gravitational collection efficiency <cgr>.
c
            vfc_1 = vf(k,i1,j1) * zmet2(1,k)
            vfc_2 = vf(k,i2,j2) * zmet2(1,k)
            cgr = e_coal * e_coll *  PI * rp**2 * abs( vfc_1 - vfc_2 )
c
c  Long's (1974) kernel that only depends on size of larger droplet
c
c           if( r_larg .le. 50.e-4 )then
c             cgr = 1.1e10 * vol(i_larg,ig_larg)**2
c           else
c             cgr = 6.33e3 * vol(i_larg,ig_larg)
c           endif
c
c  Now combine all the coagulation and collection kernels into the
c  overall kernel.
c
!            ckernel(k,i1,i2,j1,j2) = cbr + ccd + cgr
!            ckernel(k,i1,i2,j1,j2) = ZERO
!
                        ckernel(k,i1,i2,j1,j2) = cbr + cgr
!			ckernel(k,i1,i2,j1,j2) = cbr !EJL no coalesence 

!            if( itype(ienconc(j1)) .eq. I_INVOLATILE ) then

c             ---CCN-cloud particle collision --> brownian + coalescence, no charging
!              if( itype(ienconc(j2)) .eq. I_VOLATILE ) then
!               if( i1 .lt. NCCNBIN )
!     $          ckernel(k,i1,i2,j1,j2) = cbr + cgr
!              endif

c             ---For CCN-CCN collisions, use brownian coagulation only with
c                a time-dependent <ckernel> (due to charging) which will be
c                updated in {newckern}
c                Note: {newckern} assumes that the source groups are group 1

!            else  !---Target group is a cloud group

c             ---cloud particle-CCN collision --> brownian + coalescence, no charging
!              if( itype(ienconc(j2)) .eq. I_INVOLATILE ) then
!                if( i2 .lt. NCCNBIN )
!     $           ckernel(k,i1,i2,j1,j2) = cbr + cgr

c             ---cloud particle-cloud particle collision --> coalescence only
!              else
!                ckernel(k,i1,i2,j1,j2) = cgr
!              endif

!            endif

c
c  To avoid generation of large, non-physical hydrometeors by
c  coagulation, cut down ckernel for large radii
c
c           if( ( r1 .gt. 0.18 .and. r2 .gt. 10.e-4 ) .or.
c    $          ( r2 .gt. 0.18 .and. r1 .gt. 10.e-4 ) ) then
c              ckernel(k,i1,i2,j1,j2) = ckernel(k,i1,i2,j1,j2) / 1.e6
c           endif

          enddo    ! second particle bin
          enddo    ! first particle bin
         endif     ! icoag ne 0 
        enddo      ! second particle group
        enddo      ! first particle group
      enddo        ! vertical level

! E.Wolf's kluge to get the coagulation to output something other than NaNs -EJL
      open(unit=99, file='carma_ckern.txt', status='unknown')
      write(99,*) 'ckern_final'
      write(99,*) ckernel(:,:,:,1,1)  !EJL           
      close(unit=99)

c
c
c  return to caller with coagulation kernels evaluated.
c
      return
      end
