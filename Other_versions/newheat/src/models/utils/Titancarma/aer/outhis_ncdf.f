      subroutine outhis_ncdf
c
c
c  @(#) outhis.f  McKie  Jan-1997
c  This routine outputs the current model state to the history file.
c  This version outputs in netCDF format.
c  
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
c  Include constants and symbols from system netcdf header file
c   (System-dependent path to this file)
c   Following stmt is intended to be automatically edited by create_outhis_ncdf
c
      include 'netcdf.inc'
c      include 'netcdf.inc'
c      include '/ua/sys.misc/netcdf.pgi_3.4/include/netcdf.inc'
c      include '/usr/local/include/netcdf.inc'
c
c
c  Declare local variables
c
      integer id
      integer ndim
      integer nc_int
      integer nc_fp
      integer it_ncdf
      integer idim(MAXVDIMS)
      integer istart(MAXVDIMS)
      integer ncount(MAXVDIMS)
      character*(80) model_tag
      character*(30) v_units
      character*(50) v_long_name
c
c
c  Ensure values of some local variables persist between calls to this routine
c
      save nc_int
      save nc_fp
      save it_ncdf
c
c
c  Define formats
c
    1 format(/,'History output # ',i6,
     $       ' at itime: ',i6,3x,'time: ',f12.2)
c
c
c  Announce entry to this routine
c
      if( DEBUG ) write(LUNOPRT,'(/,a)') 'Enter outhis_ncdf'
c
c
c  Compute other constants
c
c
c
c  Output history header info if this is 1st call to this routine in this run
c
      if( khist .eq. 0 )then
c
c
c  Define netcdf integer type based on precision.h request symbol
c   (Allows netcdf include file to only have to be included in this routine)
c
       if( NCDF_INT .eq. NCDF_LONG )then
        nc_int = NCLONG
       else
        nc_int = NCSHORT
       endif 
c
c
c  Define netcdf floating point type based on precision.h request symbol
c   (Allows netcdf include file to only have to be included in this routine)
c
       if( NCDF_FP .eq. NCDF_DOUBLE )then
        nc_fp = NCDOUBLE
       else
        nc_fp = NCFLOAT
       endif 
c
c
c  Construct this program's id string
c
       model_tag = PROGNAM // ' (' // PROGTAG // ')'
c
c
c  Open the new netcdf file
c
       ncdf_file = nccre(hisofil, NCCLOB, ierr)
c
c
c  Declare all the netcdf variables' dimensions for variables to be in the file
c   (Do not remember all the netcdf dimension ids here)
c

       id = ncddef(ncdf_file, 'NX', NX, ierr)
       id = ncddef(ncdf_file, 'NY', NY, ierr)
       id = ncddef(ncdf_file, 'NZ', NZ, ierr)
       id = ncddef(ncdf_file, 'NXYZ', NXYZ, ierr)
       id = ncddef(ncdf_file, 'NZP1', NZ+1, ierr)

       id = ncddef(ncdf_file, 'NBIN', NBIN, ierr)
       id = ncddef(ncdf_file, 'NELEM', NELEM, ierr)
       id = ncddef(ncdf_file, 'NGROUP', NGROUP, ierr)

       id = ncddef(ncdf_file, 'NGROUPM1', NGROUPM1, ierr)

       id = ncddef(ncdf_file, 'NGAS', NGAS, ierr)

       id = ncddef(ncdf_file, 'NELPGS', NELPGS, ierr)

       id = ncddef(ncdf_file, 'NSOLUTE', NSOLUTE, ierr)

       id = ncddef(ncdf_file, 'MAXTIME', NCUNLIM, ierr)

       id = ncddef(ncdf_file, 'LEN_MODEL_TAG', len(model_tag), ierr)
       id = ncddef(ncdf_file, 'LEN_SIMTITLE', len(simtitle), ierr)
       id = ncddef(ncdf_file, 'LEN_GROUPNAME', len(groupname(1)), ierr)
       id = ncddef(ncdf_file, 'LEN_ELEMNAME', len(elemname(1)), ierr)
       id = ncddef(ncdf_file, 'LEN_GASNAME', len(gasname(1)), ierr)

c
c
c  Declare all the netcdf variables that will be output into the history file
c   (Do not try to remember all the netcdf variable ids here -- too many)
c

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'LEN_MODEL_TAG', ierr)
       id = ncvdef(ncdf_file, 'model_tag', NCCHAR, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Model name & version number'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'LEN_SIMTITLE', ierr)
       id = ncvdef(ncdf_file, 'simtitle', NCCHAR, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Simulation title'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'ibtime', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Beginning timestep number'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'ietime', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Ending timestep number'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'nhist', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of timesteps between history output'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NX', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of grid nodes in x (east-west) direction'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NY', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of grid nodes in y (north-south) direction'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NZ', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of grid nodes in the vertical'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NBIN', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of radius bins'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NELEM', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Total number of aerosol elements in all groups'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NGROUP', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of aerosol groups'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NGAS', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of gases'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 0
       id = ncvdef(ncdf_file, 'NSOLUTE', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of solutes'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'NELEM', ierr)
       id = ncvdef(ncdf_file, 'itype', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'index of aerosol type for each element'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'NELEM', ierr)
       id = ncvdef(ncdf_file, 'igelem', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Group index for each aerosol element'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'nelemg', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Number of aerosol elements for each group'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'ncore', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = '# core elements in each group'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'ienconc', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Particle # concentration element in each group'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'LEN_GROUPNAME', ierr)
       idim(2) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'groupname', NCCHAR, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Name of each aerosol group'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'LEN_ELEMNAME', ierr)
       idim(2) = ncdid(ncdf_file, 'NELEM', ierr)
       id = ncvdef(ncdf_file, 'elemname', NCCHAR, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Name of each aerosol element'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'LEN_GASNAME', ierr)
       idim(2) = ncdid(ncdf_file, 'NGAS', ierr)
       id = ncvdef(ncdf_file, 'gasname', NCCHAR, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Name of each gas'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'NBIN', ierr)
       idim(2) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'r', nc_fp, ndim, idim, ierr)
       v_units = 'cm'
       v_long_name = 'Particle radius centered at each bin'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'NBIN', ierr)
       idim(2) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'dr', nc_fp, ndim, idim, ierr)
       v_units = 'cm'
       v_long_name = 'Width of each radius bin'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'NBIN', ierr)
       idim(2) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'rmass', nc_fp, ndim, idim, ierr)
       v_units = 'g'
       v_long_name = 'Particle mass centered on each bin'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'NBIN', ierr)
       idim(2) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'dm', nc_fp, ndim, idim, ierr)
       v_units = 'g'
       v_long_name = 'width of each mass bin'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'zc', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'Altitude at middle of each vertical layer'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZP1', ierr)
       id = ncvdef(ncdf_file, 'zl', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'Altitude at boundaries of vertical layers'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'xc', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'X direction coord at center of each grid box'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'yc', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'Y direction coord at center of each grid box'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'dx', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'Grid box size in X direction'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'dy', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'Grid box size in Y direction'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'dz', nc_fp, ndim, idim, ierr)
       v_units = '(coord dep)'
       v_long_name = 'Grid box size in Z direction'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'u', nc_fp, ndim, idim, ierr)
       v_units = 'cm/s'
       v_long_name = 'X direction wind velocity at center of grid box'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'v', nc_fp, ndim, idim, ierr)
       v_units = 'cm/s'
       v_long_name = 'Y direction wind velocity at center of grid box'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NZP1', ierr)
       idim(2) = ncdid(ncdf_file, 'NBIN', ierr)
       idim(3) = ncdid(ncdf_file, 'NGROUP', ierr)
       id = ncvdef(ncdf_file, 'vf', nc_fp, ndim, idim, ierr)
       v_units = 'cm/s'
       v_long_name = 'particle fall velocity'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'itime', nc_int, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Simulation timestep counter'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'dtime', nc_fp, ndim, idim, ierr)
       v_units = 'seconds'
       v_long_name = 'Simulation time change'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 1
       idim(1) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'time', nc_fp, ndim, idim, ierr)
       v_units = 'seconds'
       v_long_name = 'Simulation time'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 6
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'NBIN', ierr)
       idim(5) = ncdid(ncdf_file, 'NELEM', ierr)
       idim(6) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'pc', nc_fp, ndim, idim, ierr)
       v_units = 'Number/cm^3'
       v_long_name = 'Aerosol particle concentration'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 5
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'NGAS', ierr)
       idim(5) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'gc', nc_fp, ndim, idim, ierr)
       v_units = 'g/cm^3'
       v_long_name = 'Gas mass concentration'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'p', nc_fp, ndim, idim, ierr)
       v_units = 'dyne/cm^2'
       v_long_name = 'Atmospheric pressure at grid box center'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       id = ncvdef(ncdf_file, 'rhoa', nc_fp, ndim, idim, ierr)
       v_units = 'g/cm^3'
       v_long_name = 'Atmospheric density at grid box center'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 4
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 't', nc_fp, ndim, idim, ierr)
       v_units = 'Deg_K'
       v_long_name = 'Temperature at center of grid box'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 4
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'ptc', nc_fp, ndim, idim, ierr)
       v_units = 'Deg_K g/cm^3'
       v_long_name = 'Potential Temperature Concentration'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 5
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'NGAS', ierr)
       idim(5) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'supsatl', nc_fp, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Supersaturation of vapor w.r.t. liquid water'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 5
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'NGAS', ierr)
       idim(5) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'supsati', nc_fp, ndim, idim, ierr)
       v_units = '(unitless)'
       v_long_name = 'Supersaturation of vapor w.r.t. ice water'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZP1', ierr)
       id = ncvdef(ncdf_file, 'w', nc_fp, ndim, idim, ierr)
       v_units = 'cm/s'
       v_long_name = 'Vertical velocity of air'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 3
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZP1', ierr)
       id = ncvdef(ncdf_file, 'dkz', nc_fp, ndim, idim, ierr)
       v_units = 'cm^2/s'
       v_long_name = 'Vertical diffusion coefficient @ vert boundary'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 5
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'NELEM', ierr)
       idim(5) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'pcmflux', nc_fp, ndim, idim, ierr)
       v_units = 'g/cm^2/s'
       v_long_name = 'Mass flux for particles'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 5
       idim(1) = ncdid(ncdf_file, 'NX', ierr)
       idim(2) = ncdid(ncdf_file, 'NY', ierr)
       idim(3) = ncdid(ncdf_file, 'NZ', ierr)
       idim(4) = ncdid(ncdf_file, 'NGAS', ierr)
       idim(5) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'gcmflux', nc_fp, ndim, idim, ierr)
       v_units = 'g/cm^2/s'
       v_long_name = 'Mass flux for gas'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'NGROUP', ierr)
       idim(2) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'inucmin', nc_int, ndim, idim, ierr)
       v_units = 'unitless'
       v_long_name = 'Minimum CN radius for nucleation'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)

       ndim = 2
       idim(1) = ncdid(ncdf_file, 'NGAS', ierr)
       idim(2) = ncdid(ncdf_file, 'MAXTIME', ierr)
       id = ncvdef(ncdf_file, 'puddle', nc_fp, ndim, idim, ierr)
       v_units = 'g/cm^2'
       v_long_name = 'Mass of liquid at surface'
       call ncaptc(ncdf_file, id, 'units', NCCHAR,
     $             len(v_units), v_units, ierr)
       call ncaptc(ncdf_file, id, 'long_name', NCCHAR,
     $            len(v_long_name), v_long_name, ierr)
c
c
c  End declarations of netcdf variables
c
       call ncendf(ncdf_file, ierr)
c
c
c
c  Output time-independent header variables to netcdf file
c

       id = ncvid(ncdf_file, 'model_tag', ierr)
       istart(1) = 1
       ncount(1) = len(model_tag)
       nlen = ncount(1)
       call ncvptc(ncdf_file, id, istart, ncount, model_tag, nlen, ierr)

       id = ncvid(ncdf_file, 'simtitle', ierr)
       istart(1) = 1
       ncount(1) = len(simtitle)
       nlen = ncount(1)
       call ncvptc(ncdf_file, id, istart, ncount, simtitle, nlen, ierr)

       id = ncvid(ncdf_file, 'ibtime', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, ibtime, ierr)

       id = ncvid(ncdf_file, 'ietime', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, ietime, ierr)

       id = ncvid(ncdf_file, 'nhist', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, nhist, ierr)

       id = ncvid(ncdf_file, 'NX', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NX, ierr)

       id = ncvid(ncdf_file, 'NY', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NY, ierr)

       id = ncvid(ncdf_file, 'NZ', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NZ, ierr)

       id = ncvid(ncdf_file, 'NBIN', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NBIN, ierr)

       id = ncvid(ncdf_file, 'NGAS', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NGAS, ierr)

       id = ncvid(ncdf_file, 'NELEM', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NELEM, ierr)

       id = ncvid(ncdf_file, 'NGROUP', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NGROUP, ierr)

       id = ncvid(ncdf_file, 'NSOLUTE', ierr)
       istart(1) = 1
       ncount(1) = 1
       call ncvpt(ncdf_file, id, istart, ncount, NSOLUTE, ierr)

       id = ncvid(ncdf_file, 'itype', ierr)
       istart(1) = 1
       ncount(1) = NELEM
       call ncvpt(ncdf_file, id, istart, ncount, itype, ierr)

       id = ncvid(ncdf_file, 'igelem', ierr)
       istart(1) = 1
       ncount(1) = NELEM
       call ncvpt(ncdf_file, id, istart, ncount, igelem, ierr)

       id = ncvid(ncdf_file, 'nelemg', ierr)
       istart(1) = 1
       ncount(1) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, nelemg, ierr)

       id = ncvid(ncdf_file, 'ncore', ierr)
       istart(1) = 1
       ncount(1) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, ncore, ierr)

       id = ncvid(ncdf_file, 'ienconc', ierr)
       istart(1) = 1
       ncount(1) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, ienconc, ierr)

       id = ncvid(ncdf_file, 'groupname', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = len(groupname(1))
       ncount(2) = NGROUP
       nlen = ncount(1) * ncount(2)
       call ncvptc(ncdf_file, id, istart, ncount, groupname, nlen, ierr)

       id = ncvid(ncdf_file, 'elemname', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = len(elemname(1))
       ncount(2) = NELEM
       nlen = ncount(1) * ncount(2)
       call ncvptc(ncdf_file, id, istart, ncount, elemname, nlen, ierr)

       id = ncvid(ncdf_file, 'gasname', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = len(gasname(1))
       ncount(2) = NGAS
       nlen = ncount(1) * ncount(2)
       call ncvptc(ncdf_file, id, istart, ncount, gasname, nlen, ierr)

       id = ncvid(ncdf_file, 'r', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = NBIN
       ncount(2) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, r, ierr)

       id = ncvid(ncdf_file, 'dr', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = NBIN
       ncount(2) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, dr, ierr)

       id = ncvid(ncdf_file, 'rmass', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = NBIN
       ncount(2) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, rmass, ierr)

       id = ncvid(ncdf_file, 'dm', ierr)
       istart(1) = 1
       istart(2) = 1
       ncount(1) = NBIN
       ncount(2) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, dm, ierr)

       id = ncvid(ncdf_file, 'zc', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, zc, ierr)

       id = ncvid(ncdf_file, 'zl', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZP1
       call ncvpt(ncdf_file, id, istart, ncount, zl, ierr)

       id = ncvid(ncdf_file, 'xc', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, xc, ierr)

       id = ncvid(ncdf_file, 'yc', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, yc, ierr)

       id = ncvid(ncdf_file, 'dx', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, dx, ierr)

       id = ncvid(ncdf_file, 'dy', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, dy, ierr)

       id = ncvid(ncdf_file, 'dz', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, dz, ierr)

       id = ncvid(ncdf_file, 'p', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, p, ierr)

       id = ncvid(ncdf_file, 'rhoa', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, rhoa, ierr)

       id = ncvid(ncdf_file, 'u', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, u, ierr)

       id = ncvid(ncdf_file, 'v', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZ
       call ncvpt(ncdf_file, id, istart, ncount, v, ierr)

       id = ncvid(ncdf_file, 'vf', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NZP1
       ncount(2) = NBIN
       ncount(3) = NGROUP
       call ncvpt(ncdf_file, id, istart, ncount, vf, ierr)

       id = ncvid(ncdf_file, 'w', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZP1
       call ncvpt(ncdf_file, id, istart, ncount, w, ierr)

       id = ncvid(ncdf_file, 'dkz', ierr)
       istart(1) = 1
       istart(2) = 1
       istart(3) = 1
       ncount(1) = NX
       ncount(2) = NY
       ncount(3) = NZP1
       call ncvpt(ncdf_file, id, istart, ncount, dkz, ierr)

c
c
c
c  Initialize netcdf timepoint counter
c   (This should be same as khist, but using a local saved var to be sure)
c
       it_ncdf = 0
c
c
c  End of time-independent header variable output
c
      endif
c
c
c  Increment netcdf timepoint counter
c
      it_ncdf = it_ncdf + 1
c
c
c  Output principal history state info to netcdf history file for current time
c

      id = ncvid(ncdf_file, 'itime', ierr)
      istart(1) = it_ncdf
      ncount(1) = 1
      call ncvpt(ncdf_file, id, istart, ncount, itime, ierr)

      id = ncvid(ncdf_file, 'dtime', ierr)
      istart(1) = it_ncdf
      ncount(1) = 1
      call ncvpt(ncdf_file, id, istart, ncount, dtime, ierr)

      id = ncvid(ncdf_file, 'time', ierr)
      istart(1) = it_ncdf
      ncount(1) = 1
      call ncvpt(ncdf_file, id, istart, ncount, time, ierr)

      id = ncvid(ncdf_file, 'pc', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = 1
      istart(5) = 1
      istart(6) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = NBIN
      ncount(5) = NELEM
      ncount(6) = 1
      call ncvpt(ncdf_file, id, istart, ncount, pc, ierr)

      id = ncvid(ncdf_file, 'gc', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = 1
      istart(5) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = NGAS
      ncount(5) = 1
      call ncvpt(ncdf_file, id, istart, ncount, gc, ierr)

      id = ncvid(ncdf_file, 't', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = 1
      call ncvpt(ncdf_file, id, istart, ncount, t, ierr)

      id = ncvid(ncdf_file, 'ptc', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = 1
      call ncvpt(ncdf_file, id, istart, ncount, ptc, ierr)

      id = ncvid(ncdf_file, 'supsatl', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = 1
      istart(5) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = NGAS
      ncount(5) = 1
      call ncvpt(ncdf_file, id, istart, ncount, supsatl, ierr)

      id = ncvid(ncdf_file, 'supsati', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = 1
      istart(5) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = NGAS
      ncount(5) = 1
      call ncvpt(ncdf_file, id, istart, ncount, supsati, ierr)

      id = ncvid(ncdf_file, 'pcmflux', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = 1
      istart(5) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = NELEM
      ncount(5) = 1
      call ncvpt(ncdf_file, id, istart, ncount, pcmflux, ierr)

      id = ncvid(ncdf_file, 'gcmflux', ierr)
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = 1
      istart(5) = it_ncdf
      ncount(1) = NX
      ncount(2) = NY
      ncount(3) = NZ
      ncount(4) = NGAS
      ncount(5) = 1
      call ncvpt(ncdf_file, id, istart, ncount, gcmflux, ierr)

      id = ncvid(ncdf_file, 'inucmin', ierr)
      istart(1) = 1
      istart(2) = it_ncdf
      ncount(1) = NGROUP
      ncount(2) = 1
      call ncvpt(ncdf_file, id, istart, ncount, inucmin, ierr)

      id = ncvid(ncdf_file, 'puddle', ierr)
      istart(1) = 1
      istart(2) = it_ncdf
      ncount(1) = NGAS
      ncount(2) = 1
      call ncvpt(ncdf_file, id, istart, ncount, puddle, ierr)
c
c  Synchronize netcdf disk output file with netcdf disk
c
      call ncsnc(ncdf_file, ierr)

c
c
c  Return to caller with history output complete
c
      return
      end