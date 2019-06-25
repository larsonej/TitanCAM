#  @(#) Makefile  McKie  Sep-1999
#  This is the pseudo make file for manipulating the dummy netcdf code.
#
#  All real work is performed by the Makefile.real make file, which
#  is intended to be called by the top Makefile, so that critical
#  symbol values defined in the top Makefile can be passed.
#
#  This Makefile therefore calls the top Makefile for each of its
#  targets, which then calls the Makefile.work in this directory
#  to do the work.
#
#  Modified to use Makefile and Makefile.work and top level Makefile
#  control, Sep-1999, for carma v1.18.  -bm
#
#
#  Begin target list
#
#  Targets:
#
#    default     Same as library
#    library     Ensure include files are consistent, compile, and build library 
#    compile     Ensure each latest source module is compiled
#    clean       Remove all compiler object files
#    tar         Build distribution tar for dummy netcdf only code
#    help        List available useful targets
#
#  End target list
#

default:	;
		( cd .. ; make -f Makefile default_dumnetcdf )

library:	;
		( cd .. ; make -f Makefile library_dumnetcdf )

compile:	;
		( cd .. ; make -f Makefile compile_dumnetcdf )

clean:		;
		( cd .. ; make -f Makefile clean_dumnetcdf )

tar:		;
		( cd .. ; make -f Makefile tar_dumnetcdf )

help:		;
		@ sed -n '/^#  Begin target list/,/^#  End target list/p'  Makefile
