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
      include 'globals.h'
c
c
c  Do model initializations
c
      call init
c
c  temperature oscillation initialization (for testing purposes only)
c
      tzero = t(1,1,1)
      tosc_amp = .1
      tosc_per = 30.
c
c  Do the requested number of timesteps for this run
c
      do i=ibtime+1,ietime

       call step
c
c  temperature oscillation (for testing purposes only)
c
       t(1,1,1) = tzero - tosc_amp*sin(2.*PI*time/tosc_per)

      enddo
c
c
c  Do normal program terminations
c
      call quit
      end
