#  @(#) Makefile  McKie  Sep-1999
#  This is the pseudo make file for manipulating the radiation submodel code.
#
#  All real work is performed by the Makefile.work make file, which
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
#    tar         Build distribution tar for radiation only code
#    help        List possible useful targets
#
#  End target list
#

default:	;
		( cd .. ; make -f Makefile default_rad )

library:	;
		( cd .. ; make -f Makefile library_rad )

compile:		;
		( cd .. ; make -f Makefile compile_rad )

clean:		;
		( cd .. ; make -f Makefile clean_rad )

tar:		;
		( cd .. ; make -f Makefile tar_rad )

help:		;
		@ sed -n '/^#  Begin target list/,/^#  End target list/p'  Makefile
