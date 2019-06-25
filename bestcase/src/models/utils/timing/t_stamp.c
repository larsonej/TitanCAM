/*
** $Id: t_stamp.c 17 2006-12-11 21:50:24Z hpc $
*/

#include <sys/time.h>     /* gettimeofday */ 
#include <sys/times.h>    /* times */
#include <unistd.h>       /* gettimeofday */

#include <gpt.h>

/*
** t_stamp: Compute timestamp of usr, sys, and wallclock time (seconds)
**
** Output arguments:
**   wall: wallclock
**   usr:  user time
**   sys:  system time
**
** Return value: 0 (success) or -1 (failure)
*/

int t_stamp (double *wall, double *usr, double *sys)
{
  struct timeval tp;         /* argument to gettimeofday */
  struct tms buf;            /* argument to times */

  *usr = 0;
  *sys = 0;

  if (times (&buf) == -1)
    return t_error ("t_stamp: times() failed. Timing bogus\n");

  *usr = buf.tms_utime / (double) ticks_per_sec;
  *sys = buf.tms_stime / (double) ticks_per_sec;

  gettimeofday (&tp, NULL);
  *wall = tp.tv_sec + 1.e-6*tp.tv_usec;

  return 0;
}
