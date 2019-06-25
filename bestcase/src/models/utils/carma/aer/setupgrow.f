       subroutine setupgrow
c
c
c  @(#) setupgrow.f  Ackerman  Dec-1995
c  This routine defines time-independent parameters used to calculate
c  condensational growth/evaporation.
c
c  The parameters defined for each gas are 
c
c    gwtmol:   molecular weight [g/mol]
c    diffus:   diffusivity      [cm^2/s]
c    rlhe  :   latent heat of evaporation [cm^2/s^2]
c    rlhm  :   latent heat of melting [cm^2/s^2]
c
c  Time-independent parameters that depend on particle radius are
c  defined in setupgkern.f.
c
c  This routine requires that vertical profiles of temperature <T>,
c  and pressure <p> are defined (i.e., initatm.f must be called before this).
c  The vertical profile with ix = iy = 1 is used.
c
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
    1 format(a,':  ',12i6)
    2 format(a,':  ',i6)
    3 format(/' id  gwtmol   gasname',(/,i3,3x,f5.1,3x,a))
    5 format(/,'Particle growth mapping arrays (setupgrow):')
c
c-------------------------------------------------------------------------------
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter setupgrow'
c
c-----Check that values are valid------------------------------------------
c
      do ielem = 1, NELEM

        if( igrowgas(ielem) .gt. NGAS )then
          write(LUNOPRT,'(/,a)') 'component of igrowgas > NGAS'
          call endcarma
        endif

      enddo
c
c--------------------------------------------------------------------------
c
c  
c  Define parameters with weak time-dependence to be used in
c  growth equation.
c
      do k = 1, NZ
c
c
c   Diffusivity of water vapor in air from Pruppacher & Klett (eq. 13-3);
c   units are [cm^2/s].
c
        diffus(k,1) =  0.211 * ( 1.01325e+6 / p(1,1,k) )
     $              * ( t(1,1,k) / 273.15 )**1.94
c
c
c   Latent heat of evaporation for water from Stull; units are [cm^2/s^2]
c
        rlhe(k,1) = ( 2.5 - .00239*( t(1,1,k) - 273.16 ) )*1.e10
c
c
c   Latent heat of ice melting from Pruppacher & Klett (eq. 4-85b);
c   units are [cm^2/s^2]
c
        rlhm(k,1) = ( 79.7 + 0.485*(t(1,1,k)-273.16) - 2.5e-3*
     1               ( (t(1,1,k) - 273.16)**2 ) )*4.186e7

      enddo
c
c
c--------------------------------------------------------------------------
c
c
c  Report some initialization values
c
      write(LUNOPRT,5)
      write(LUNOPRT,2) 'NGAS    ',NGAS
      write(LUNOPRT,1) 'igrowgas',(igrowgas(i),i=1,NELEM)
      write(LUNOPRT,3) (i,gwtmol(i),gasname(i),i=1,NGAS)
c
c
c  Return to caller with particle growth mapping arrays and time-dependent
c  parameters initialized.
c
      return
      end
