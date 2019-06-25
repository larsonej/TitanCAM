      subroutine initres
c
c
c  @(#) initres.f  McKie  Oct-1995
c  This routine initializes the model from a restart that
c  was created in a previous run.  Common blocks are
c  returned to the state they were in when the requested
c  restart time values were output in the previous run.
c
c  Argument list input:
c    None.
c
c  Argument list output:
c    None.
c
c
c  Include global constants and variables
c
      include 'globaer.h'
c
c
c  Define formats
c
    1 format(/,'Attempting to restart from file: ',a)
    2 format('Looking for restart time index: ',i6)
    3 format('Unexpected common block value for ',a,': ',i6)
    4 format('Unexpected common block value for ',a,': ',a)
    5 format('Input restart file time:  itime: ',i6,3x,'time: ',f12.2)
    7 format('Requested restart time index ',i6,' found.')
    8 format('Error--Cant find requested restart time ',i6,
     $       ' on input restart file.')
    9 format(/,'Doing a restart initialization from a prev run')
   10 format('Aborting.')
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter initres'
c
c
c  Report that a restart is being done
c
      write(LUNOPRT,9)
c
c
c  Report input restart file name
c
      call dblank(resifil, ns)
      write(LUNOPRT,1) resifil(1:ns)
c
c
c  Report requested restart timestep index
c
      write(LUNOPRT,2) ibtime
c
c
c  Open input restart file
c
      open(unit=LUNIRES,file=resifil,status='old',
     $     form='unformatted')
c
c
c  Input restart records for next time on restart file
c
 2100 continue
c
c
c
c
c  Input contents of common /aer1/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    dom_llx, dom_lly, dom_urx, dom_ury, rlon0, rlat0, rlat1,
     $    rlat2, hemisph, zl, zc, zlold, zcold, xc, xl, xu, yc, yl, yu,
     $    dx, dy, dz, xmet, ymet, zmet, zmetl, rlon, rlat, igridv,
     $    igridh, iaer1
       if( iaer1 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer1', iaer1
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer1s/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    gridname, caer1s
       if( caer1s .ne. CSAFETY )then
        write(LUNOPRT,4) 'caer1s', caer1s
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer2/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    time, dtime, dtmin, dtmax, dpctol, dgstol, conmax, time_nuc,
     $    period_nuc, maxsubsteps, minsubsteps, dtime_save, sec_mom,
     $    rlh_nuc, itime, ix, iy, iz, ixy, ixyz, ntsubsteps, igelem,
     $    itype, icomp, nelemg, ncore, ishape, ienconc, icoag,
     $    icoagelem, ifall, icoagop, icollec, itbnd, ibbnd, itbnd_pc,
     $    ibbnd_pc, itbnd_gc, ibbnd_gc, itbnd_ptc, ibbnd_ptc, ihoradv,
     $    do_coag, do_grow, do_thermo, do_vtran, do_ew, do_ns,
     $    do_varstep, do_step, do_ccoef, do_error, do_netcdf, do_parcel,
     $    if_nuc, if_nuc_lh, ncdf_file, imomelem, inucproc, igrowgas,
     $    inucgas, nnuc2elem, inuc2elem, ievp2elem, is_grp_ice,
     $    is_grp_mixed, inuc2bin, isolelem, ievp2bin, icorelem,
     $    nnucelem, nnucbin, inucelem, inucbin, iaer2
       if( iaer2 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer2', iaer2
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer2s/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    simtitle, elemname, groupname, gasname, solname, caer2s
       if( caer2s .ne. CSAFETY )then
        write(LUNOPRT,4) 'caer2s', caer2s
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer3/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    rmin, rmassmin, rmrat, r, rmass, vol, dr, dm, dv, rmassup,
     $    rup, rlow, diffmass, rhop, rhoelem, eshape, iaer3
       if( iaer3 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer3', iaer3
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer4/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    pc, gc, ptc, iaer4
       if( iaer4 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer4', iaer4
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer5/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    pcl, gcl, ptcl, d_pc, d_gc, d_ptc, pcmax, cvert, divcor, chor,
     $    dhor, pconmax, coaglg, coagpe, rnuclg, rnucpe, growlg, growpe,
     $    evaplg, evappe, gasprod, rlheat, vertdifd, vertdifu, ftopgas,
     $    fbotgas, ftoppart, fbotpart, ftop, fbot, cmf, totevap,
     $    inucmin, inucstep, rprod, ppd, pls, iaer5
       if( iaer5 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer5', iaer5
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer6/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    ck0, grav_e_coll0, ckernel, pkernel, volx, ilow, jlow, iup,
     $    jup, npairl, npairu, cbr_term0, cbr_term1, cbr_term2, cgr,
     $    pkern0, iaer6
       if( iaer6 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer6', iaer6
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer7/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    iglow, jglow, igup, jgup, iaer7
       if( iaer7 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer7', iaer7
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer8/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    bpm, vf, vtrans, re, vf_const, vertadvu, vertadvd, htrans,
     $    hdiff, ca, cb, cd, ce, cf, cg, iaer8
       if( iaer8 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer8', iaer8
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer9/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    zbot, p, rhoa, p_surf, p_top, t, told, pold, rhoaold, rmu,
     $    thcond, w, zmetold, u, v, t_surf, dkz, dkx, dky, iaer9
       if( iaer9 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer9', iaer9
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer10/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    gwtmol, diffus, rlhe, rlhm, pvapl, pvapi, surfctwa, surfctiw,
     $    surfctia, akelvin, akelvini, ft, gro, gro1, gro2, gvrat,
     $    supsatl, supsati, supsatlold, supsatiold, scrit, sol_ions,
     $    solwtmol, rhosol, iaer10
       if( iaer10 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer10', iaer10
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /aer11/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    adelf, bdelf, prenuc, rmiv, iaer11
       if( iaer11 .ne. ISAFETY )then
        write(LUNOPRT,3) 'iaer11', iaer11
        write(LUNOPRT,10)
        stop
       endif
c
c
c  Input contents of common /rad1/ & check for safety marker value
c
       read(LUNIRES,end=8100,err=8100)
     $    u0, u0_fixed, rad_start, zsin, zcos, wave, qrad, radheat,
     $    alb_tomi, alb_toai, alb_toa, opd, fsl_up, fsl_dn, fir_up,
     $    fir_dn, do_rad, do_solar, do_ir, nrad, prad, isolar_zen, irad1
       if( irad1 .ne. ISAFETY )then
        write(LUNOPRT,3) 'irad1', irad1
        write(LUNOPRT,10)
        stop
       endif
c
c
c
c  If current restart file time is not the requested one, go try next time
c
      if( itime .ne. ibtime ) goto 2100
c
c
c  Close input restart file
c
      close(unit=LUNIRES)
c
c
c  Report that requested restart info has been input
c
      call dblank(resifil, ns)
      write(LUNOPRT,7) itime
c
c
c  Return to caller with restart initialization complete
c
      return
c
c
c  Report that requested restart time was not found on input file & quit
c
 8100 continue
      write(LUNOPRT,8) ibtime
      write(LUNOPRT,10)
      stop
c
      end
