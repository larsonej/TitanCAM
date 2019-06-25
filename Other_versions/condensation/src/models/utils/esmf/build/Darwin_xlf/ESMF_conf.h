#ifdef ESMC_RCS_HEADER
"$Id: ESMF_conf.h 62 2008-04-23 22:59:18Z cam_titan $"
"Defines the configuration for this machine"
#endif

#if !defined(INCLUDED_CONF_H)
#define INCLUDED_CONF_H

#define PARCH_mac_osx
#define ESMF_ARCH_NAME "mac_osx"

#define ESMC_HAVE_LIMITS_H
#define ESMC_HAVE_PWD_H 
#define ESMC_HAVE_MALLOC_H 
#define ESMC_HAVE_STRING_H 
#define ESMC_HAVE_GETDOMAINNAME
#define ESMC_HAVE_DRAND48 
#define ESMC_HAVE_UNAME 
#define ESMC_HAVE_UNISTD_H 
#define ESMC_HAVE_SYS_TIME_H 
#define ESMC_HAVE_STDLIB_H

#define ESMC_HAVE_FORTRAN_UNDERSCORE
#define ESMC_HAVE_FORTRAN_UNDERSCORE_UNDERSCORE

#define ESMC_POINTER_SIZE 8
#undef ESMC_HAVE_OMP_THREADS

#define ESMC_HAVE_MPI 1

#define ESMC_HAVE_READLINK
#define ESMC_HAVE_MEMMOVE
#define ESMC_HAVE_TEMPLATED_COMPLEX

#define ESMC_HAVE_DOUBLE_ALIGN_MALLOC
#define ESMC_HAVE_MEMALIGN
#define ESMC_HAVE_SYS_RESOURCE_H
#define ESMC_SIZEOF_VOIDP 4
#define ESMC_SIZEOF_INT 4
#define ESMC_SIZEOF_DOUBLE 8

#define ESMC_HAVE_SYS_UTSNAME_H
#endif
