#include <misc.h>

!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module tmp_qrs_mod

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver

  implicit none
  public          ! Make default type public 

  real(r8) ::  tmp_qrs(pcols,pver)
 
end module tmp_qrs_mod 
