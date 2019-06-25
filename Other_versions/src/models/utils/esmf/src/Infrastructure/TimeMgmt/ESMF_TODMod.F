! $Id: ESMF_TODMod.F 62 2008-04-23 22:59:18Z cam_titan $
	module ESMF_TODMod
!===============================================================================
!BOP
! !MODULE: ESMF_TODMod
! 
! !USES:
!
! !PUBLIC TYPES:
	implicit none

	type ESMF_TOD 
          private
          sequence
          integer(8) type
          integer(8) sec
          integer(8) msec
	end type ESMF_TOD 
!
! !DESCRIPTION:
! Describes and contains a time of day.
!EOP
!===============================================================================

	end module
