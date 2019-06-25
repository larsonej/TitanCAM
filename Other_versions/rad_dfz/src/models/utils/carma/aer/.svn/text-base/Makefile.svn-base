#  @(#) Makefile  McKie  Sep-1999
#  This is the pseudo make file for manipulating the aerosol submodel code.
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
#    default     Call top level make to generate executable
#    library     Ensure include files are consistent, compile, and build library 
#    compile     Ensure each latest source module is compiled
#    clean       Remove all compiler object files
#    do_chekres  Build & run program to check restart i/o source code against commons
#    tar         Build distribution tar for aerosol only code
#    help        List available useful targets
#
#  End target list
#

default:	;
		( cd .. ; make -f Makefile default )

library:	;
		( cd .. ; make -f Makefile library_aer )

compile:		;
		( cd .. ; make -f Makefile compile_aer )

clean:		;
		( cd .. ; make -f Makefile clean_aer )

tar:		;
		( cd .. ; make -f Makefile tar_aer )

do_chekres:	;
		( cd .. ; make -f Makefile do_chekres_aer )

help:		;
		@ sed -n '/^#  Begin target list/,/^#  End target list/p'  Makefile
