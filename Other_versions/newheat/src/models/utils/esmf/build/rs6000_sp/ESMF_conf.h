#ifdef ESMC_RCS_HEADER
"$Id: ESMF_conf.h 17 2006-12-11 21:50:24Z hpc $"
"Defines the configuration for this machine"
#endif

#if !defined(INCLUDED_CONF_H)
#define INCLUDED_CONF_H

#define PARCH rs6000_sp
#define ESMF_ARCH_NAME "rs6000_sp"
#define ESMC_USE_READ_REAL_TIME

#define ESMC_HAVE_LIMITS_H
#define ESMC_HAVE_STROPTS_H 
#define ESMC_HAVE_SEARCH_H 
#define ESMC_HAVE_PWD_H 
#define ESMC_HAVE_STDLIB_H
#define ESMC_HAVE_STRING_H 
#define ESMC_HAVE_STRINGS_H 
#define ESMC_HAVE_MALLOC_H 
#define _POSIX_SOURCE 
#define ESMC_HAVE_DRAND48  
#define ESMC_HAVE_GETDOMAINNAME 
#if !defined(_XOPEN_SOURCE)
#define _XOPEN_SOURCE 
#endif
#define ESMC_HAVE_UNISTD_H 
#define ESMC_HAVE_SYS_TIME_H 
#define ESMC_HAVE_UNAME 
#if !defined(_XOPEN_SOURCE_EXTENDED)
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#define _ALL_SOURCE   
#define ESMC_HAVE_BROKEN_REQUEST_FREE 
#define ESMC_HAVE_STRINGS_H
#define ESMC_HAVE_DOUBLE_ALIGN_MALLOC

#undef ESMC_HAVE_MPI

#define ESMC_SUBSTITUTE_CTRL_CHARS 1

#define ESMC_POINTER_SIZE 8

#define ESMC_HAVE_XLF90

#define ESMC_PREFER_BZERO

#define ESMC_HAVE_READLINK
#define ESMC_HAVE_MEMMOVE

#define ESMC_HAVE_PRAGMA_DISJOINT

#define ESMC_USE_DBX_DEBUGGER
#define ESMC_HAVE_SYS_RESOURCE_H
#define ESMC_SIZEOF_VOIDP 4
#define ESMC_SIZEOF_INT 4
#define ESMC_SIZEOF_DOUBLE 8

#define ESMC_WORDS_BIGENDIAN 1
#define ESMC_NEED_SOCKET_PROTO
#define ESMC_HAVE_ACCEPT_SIZE_T

#define ESMC_HAVE_SYS_UTSNAME_H

#define ESMC_HAVE_SLEEP_RETURNS_EARLY
#define ESMC_USE_KBYTES_FOR_SIZE

#define ESMC_USE_A_FOR_DEBUGGER

#endif
