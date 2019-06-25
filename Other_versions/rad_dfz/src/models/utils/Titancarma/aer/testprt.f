       subroutine testprt
c
c
c  @(#) testprt.f  Barth  Feb-2003
c  This routine controls which write statements will be executed when 
c  debugging
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
c
c  Altitude for outputs
      iwa = 21
      if(ext .eq. 'ge') iwa = 26
      if(ext .eq. 'gh') iwa = 8
      if(ext .eq. 'gi') iwa = 8

c  Bin for outputs
      iwb = 20

c  {csolve} 
      do_write(12) = .false.   !f>1 after coag

c  {growevapl}
      do_write(1) = .false.   !growcore mass fraction
      do_write(2) = .false.   !growth matrix terms
      do_write(11) = .false.  !Courant # > 1
 
c  {psolve}
      do_write(3) = .true.   !mass conc of cloud and growcore
      do_write(10) = .false.   !cloud g/cm3 for each bin
    
c  {gasexchange}
      do_write(4) = .false.   !gas production terms
      do_write(14) = .true. !gprod_grow terms
 
c  {downgevapply}
      do_write(8) = .true.   !Nucl/Evap production terms

c  {isotope_exc}
      do_write(5) = .true.   !gas exchange terms
 
c  {gsolve}
      do_write(6) = .true.   !new gc

c  {varstep} 
      do_write(7) = .true.   !negative gas concentrations

c  Print pc/gc after several subroutines
      do_write(9) = .false.

c  Print fractions (just set in newstate for now)
      do_write(15) = .true.

      return
      end
