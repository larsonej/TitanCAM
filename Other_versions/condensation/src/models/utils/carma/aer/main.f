      program main
c
c
c  @(#) main.f  McKie  Oct-1995
c  This is the main program for the Ames Aerosol Model.
c  If this model is being used as a submodel from another
c  program, the things done in this program will be
c  done in the master model's main program.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Do model initializations
c
      call init
c
c
c  Take timesteps until requested number of timesteps taken
c  or simulation time exceeds maximum, whichever comes first.
c
 1100 if( (itime.le.ietime) .and. (time.lt.(endtime-.5*dtime)) )then
       call step
       goto 1100
      endif
c
c
c  Do normal program terminations
c
      call quit
      end
