/*
 * $Id: scam.c 17 2006-12-11 21:50:24Z hpc $
 *  main() routine for fifo based client-server version
 */

#include <rpc/rpc.h>
#include <sys/time.h>

#include "fortran.h"
#include "scam_fifo.h"
#ifndef HP
#include "scam_rpc.h"
#endif

#ifndef lint
static char rcsid[] = "$Id: scam.c 17 2006-12-11 21:50:24Z hpc $";
#endif 

#ifndef DEBUG_INIT
int
main( int argc, char* argv[] )
{
    server_init( argc, argv );
    return server_run();
}
    
#endif









