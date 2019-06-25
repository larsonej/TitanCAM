      subroutine outres
c
c
c  @(#) outres.f  McKie  Oct-1995
c  This routine outputs the current critical global model
c  information to the output restart file.
c  This routine is created automatically from a template file.
c  So edits should be made to the template file, otherwise they
c  could be lost due to recreation.
c  (Template file is usually outres.template)
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
    1 format('Restart info written at itime = ',i6,
     $       3x,'time=',f12.4,4x,'to file ',a)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter outres'
c
c
c
c
c  Load safety value & output contents of common /aer1/
c
       iaer1 = ISAFETY
       write(LUNORES)
     $    dom_llx, dom_lly, dom_urx, dom_ury, rlon0, rlat0, rlat1,
     $    rlat2, hemisph, zl, zc, zlold, zcold, xc, xl, xu, yc, yl, yu,
     $    dx, dy, dz, xmet, ymet, zmet, zmetl, rlon, rlat, igridv,
     $    igridh, iaer1
c
c
c  Load safety value & output contents of common /aer1s/
c
       caer1s = CSAFETY
       write(LUNORES)
     $    gridname, caer1s
c
c
c  Load safety value & output contents of common /aer2/
c
       iaer2 = ISAFETY
       write(LUNORES)
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
c
c
c  Load safety value & output contents of common /aer2s/
c
       caer2s = CSAFETY
       write(LUNORES)
     $    simtitle, elemname, groupname, gasname, solname, caer2s
c
c
c  Load safety value & output contents of common /aer3/
c
       iaer3 = ISAFETY
       write(LUNORES)
     $    rmin, rmassmin, rmrat, r, rmass, vol, dr, dm, dv, rmassup,
     $    rup, rlow, diffmass, rhop, rhoelem, eshape, iaer3
c
c
c  Load safety value & output contents of common /aer4/
c
       iaer4 = ISAFETY
       write(LUNORES)
     $    pc, gc, ptc, iaer4
c
c
c  Load safety value & output contents of common /aer5/
c
       iaer5 = ISAFETY
       write(LUNORES)
     $    pcl, gcl, ptcl, d_pc, d_gc, d_ptc, pcmax, cvert, divcor, chor,
     $    dhor, pconmax, coaglg, coagpe, rnuclg, rnucpe, growlg, growpe,
     $    evaplg, evappe, gasprod, rlheat, vertdifd, vertdifu, ftopgas,
     $    fbotgas, ftoppart, fbotpart, ftop, fbot, cmf, totevap,
     $    inucmin, inucstep, rprod, ppd, pls, iaer5
c
c
c  Load safety value & output contents of common /aer6/
c
       iaer6 = ISAFETY
       write(LUNORES)
     $    ck0, grav_e_coll0, ckernel, pkernel, volx, ilow, jlow, iup,
     $    jup, npairl, npairu, cbr_term0, cbr_term1, cbr_term2, cgr,
     $    pkern0, iaer6
c
c
c  Load safety value & output contents of common /aer7/
c
       iaer7 = ISAFETY
       write(LUNORES)
     $    iglow, jglow, igup, jgup, iaer7
c
c
c  Load safety value & output contents of common /aer8/
c
       iaer8 = ISAFETY
       write(LUNORES)
     $    bpm, vf, vtrans, re, vf_const, vertadvu, vertadvd, htrans,
     $    hdiff, ca, cb, cd, ce, cf, cg, iaer8
c
c
c  Load safety value & output contents of common /aer9/
c
       iaer9 = ISAFETY
       write(LUNORES)
     $    zbot, p, rhoa, p_surf, p_top, t, told, pold, rhoaold, rmu,
     $    thcond, w, zmetold, u, v, t_surf, dkz, dkx, dky, iaer9
c
c
c  Load safety value & output contents of common /aer10/
c
       iaer10 = ISAFETY
       write(LUNORES)
     $    gwtmol, diffus, rlhe, rlhm, pvapl, pvapi, surfctwa, surfctiw,
     $    surfctia, akelvin, akelvini, ft, gro, gro1, gro2, gvrat,
     $    supsatl, supsati, supsatlold, supsatiold, scrit, sol_ions,
     $    solwtmol, rhosol, iaer10
c
c
c  Load safety value & output contents of common /aer11/
c
       iaer11 = ISAFETY
       write(LUNORES)
     $    adelf, bdelf, prenuc, rmiv, iaer11
c
c
c  Load safety value & output contents of common /rad1/
c
       irad1 = ISAFETY
       write(LUNORES)
     $    u0, u0_fixed, rad_start, zsin, zcos, wave, qrad, radheat,
     $    alb_tomi, alb_toai, alb_toa, opd, fsl_up, fsl_dn, fir_up,
     $    fir_dn, do_rad, do_solar, do_ir, nrad, prad, isolar_zen, irad1
c
c
c
c  Report that restart info has been output
c
      call prtsep
      call dblank(resofil, ns)
      write(LUNOPRT,1) itime, time, resofil(1:ns)
c
c
c  Return to caller with model info output to restart file
c
      return
      end
