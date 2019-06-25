c  @(#) project.f  McKie  Jun=1997
c  This collection of routines handles geographical projections
c  and calcs associated with those projections.
c
c  For each projection, there are usually 3 routines:
c    Compute fixed parameters of the projection.
c    Compute the forward projection from (lon,lat) to projection plane (u,v).
c    Compute the inverse projection from projection plane (u,v) to (lon,lat).
c
c  Projections currently handled:
c    (lc)  Lambert Conformal
c    (ps)  Polar Stereographic
c    (me)  Mercator
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      function sfirdis(rlon1,rlat1,rlon2,rlat2)
c
c
c  @(#) sfirdis  McKie  Jul-1988
c  This routine computes the great circle distance between 2 given
c  points on the surface of a unit sphere.
c
c  Input:
c          rlon1,rlat1 = longitude,latitude (degrees) of point 1
c          rlon2,rlat2 = longitude,latitude (degrees) of point 2
c
c  Returns:  Shortest great circle distance between pts 1 & 2 in radians
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram arg(s)
c
c     real rlon1,rlat1
c     real rlon2,rlat2
c
c
c  Define degree to radian factor
c
      parameter(DEG2RAD=.0174532925199433)
c
c
c  Compute radian form of input angles
c
      a1 = DEG2RAD*rlat1
      b1 = DEG2RAD*rlon1
      a2 = DEG2RAD*rlat2
      b2 = DEG2RAD*rlon2
c
c
c  Compute cos of angle on great circle separating pts 1 & 2
c   (dot prod of the 2 vectors from sphere centre to pts 1 & 2)
c
      cosine = cos(a1)*cos(a2) * ( cos(b1)*cos(b2) + sin(b1)*sin(b2) )
     $              + sin(a1)*sin(a2)
c
c
c  Return to caller with distance measure in value of this function
c
      sfirdis = acos(cosine)
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
       subroutine parmlc(rlat1,rlat2, sgn,sgn2,cone)
c
c
c  @(#) parmlc.f  McKie Aug-1989
c  This routine computes parameters associated with the lambert
c  conformal conic projection with 2 standard parallels.
c
c  Input:
c    rlat1,rlat2 = 2 standard parallels
c
c  Output:
c            sgn = hemisphere constant: +1.=southern, -1.=northern
c           sgn2 = .5 * sgn
c           cone = cone constant
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlat1,rlat2
c     real sgn,sgn2
c     real cone
c
c
c  Define local symbolic constants
c
      parameter(DEG2RAD=.0174532925199433)
c
c
c  Compute hemisphere constants
c
      if( .5 * ( rlat1 + rlat2 ) .ge. 0. )then
       sgn = -1.
      else
       sgn = 1.
      endif
      sgn2 = .5 * sgn
c
c
c  Compute cone constant
c
      a1 = log( cos( DEG2RAD*rlat1 ) ) 
      a2 = log( cos( DEG2RAD*rlat2 ) ) 
      b1 = log( tan( DEG2RAD*(45.+sgn2*rlat1) ) )
      b2 = log( tan( DEG2RAD*(45.+sgn2*rlat2) ) )
      cone = ( a1 - a2 ) / ( b1 - b2 )
c
c
c  Return to caller with projection parameters
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
       subroutine parmme(rlat0, factme)
c
c
c  @(#) parmme.f  McKie  Mar-1992
c  This routine computes a parameter associated with the 
c  Mercator projection.  (This routine is usually called
c  once to precompute the constant factme used in succeeding
c  calls to projme & invpme.)
c
c  Input:
c          rlat0 = latitude at which projection is true (proj plane here)
c
c  Output:
c         factme = projection factor
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlat0
c     real factme
c
c
c  Define local symbolic constants
c
      parameter(DEG2RAD=.0174532925199433)
      parameter(POLE=89.99)
c
c
c  Ensure that factor doesn't lead to trivial projection into (0,0)
c
      angle = DEG2RAD * max( -POLE, min( POLE, rlat0 ) )
c
c
c  Compute projection factor
c
      factme = cos( angle )
c
c
c  Return to caller with projection parameter
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
       subroutine parmps(rlat0,sgn, factps)
c
c
c  @(#) parmps.f  McKie  Sep-1991
c  This routine computes a parameter associated with the 
c  polar stereographic projection.  (This routine is
c  usually called once to precompute the constant used
c  in succeeding calls to projps & invpps.)
c
c  Input:
c          rlat0 = latitude at which projection is true (proj plane here)
c            sgn = hemisphere constant: +1.=southern, -1.=northern
c
c  Output:
c         factps = projection factor
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlat0
c     real sgn
c     real factps
c
c
c  Define local symbolic constants
c
      parameter(DEG2RAD=.0174532925199433)
c
c
c  Ensure sgn is reasonable
c
      if( (abs(sgn+1.) .gt. .0001) .and. (abs(sgn-1.) .gt. .0001))then
       write(*,*) 'Error--parmps: sgn=',sgn,' not + or - 1.'
       stop
      endif
c
c
c  Compute projection factor
c
      factps = abs( sin(DEG2RAD*rlat0) - sgn )
c
c
c  Check if factor is reasonable
c
      if( factps .lt. .00001 )then
       write(*,*) 'Error--parmps: projection plane on projection pole'
       write(*,*) '               rlat0,sgn=',rlat0,sgn
       stop
      endif
c
c
c  Return to caller with projection parameter
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
       subroutine projlc(cone,rlon0,sgn,sgn2,n, rlond,rlatd, u,v)
c
c
c  @(#) projlc.f  McKie Jul-1989
c  This routine performs the Lambert conformal projection with 2 standard
c  parallels from longitude & latitude coordinates to u,v space.
c
c  Input:
c          cone               = cone constant
c          rlon0              = central longitude
c          sgn                = -1. for northern hemisphere, +1. for southern
c          sgn2               = sgn / 2.
c          n                  = # coord pairs to project
c          rlatd(i),rlond(i)  = lon,lat coord (in degrees)
c
c  Output:
c          u(i),v(i)        = projected coord 
c
c  Note:  rlatd(i) & u(i) can safely be same array, as can be rlond(i) & v(i).
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real cone
c     real rlon0
c     real sgn
c     real sgn2
c     integer n
      dimension rlatd(n),rlond(n)
      dimension u(n),v(n)
c
c
c  Define degrees to radians factor symbolic constant
c
      parameter(DEG2RAD=.0174532925199433)
      parameter(POLE=89.99)
c
c
c  Compute a constant factor
c
      factor = DEG2RAD * cone
c
c
c  Compute projected coordinates for each input longitude & latitude
c
      do 3100 i=1,n
       y = max( -POLE, min( POLE, rlatd(i) ) )
       r = ( tan( DEG2RAD * ( 45. + sgn2*y ) ) ) ** cone
       azimuth = factor * ( rlond(i) - rlon0 )
       u(i) = r * sin(azimuth)
       v(i) = sgn * r * cos(azimuth)
 3100 continue
c
c
c  Return to caller with u(i),v(i)
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine projme(rlon0,factme, n, rlond,rlatd, u,v)
c
c
c  @(#) projme.f  McKie  Mar-1992
c  This routine performs the Mercator projection 
c  from longitude, latitude coordinates to u,v space.
c
c  Input:
c          rlon0              = central longitude
c          factme             = Mercator factor, cos(rlat0)
c          n                  = # coord pairs to project
c          rlatd(i),rlond(i)  = lon,lat coord (in degrees)
c
c  Output:
c          u(i),v(i)        = projected coord 
c
c  Note:  rlond(i) & u(i) can safely be same array, as can be rlatd(i) & v(i).
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlon0
c     real factme
c     integer n
      dimension rlatd(n),rlond(n)
      dimension u(n),v(n)
c
c
c  Define degrees to radians factor symbolic constant
c
      parameter(DEG2RAD=.0174532925199433)
      parameter(POLE=89.99)
c
c
c  Compute projected coordinates for each input longitude & latitude
c
      do i=1,n
       phi = DEG2RAD * max( -POLE, min( POLE, rlatd(i) ) )
       u(i) = factme * DEG2RAD * ( rlond(i) - rlon0 )
       v(i) = factme * log( ( 1. + sin(phi) ) / cos(phi) ) 
      enddo
c
c
c  Return to caller with u(i),v(i)
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine projps(rlon0,sgn,factps, n, rlond,rlatd, u,v)
c
c
c  @(#) projps.f  McKie Sep-1991
c  This routine performs the polar stereographic projection 
c  from longitude, latitude coordinates to u,v space.
c
c  Input:
c          rlon0              = central longitude
c          sgn                = -1. for northern hemisphere, +1. for southern
c          factps             = polar stereographic factor, abs(sin(phi0) - sgn)
c          n                  = # coord pairs to project
c          rlatd(i),rlond(i)  = lon,lat coord (in degrees)
c
c  Output:
c          u(i),v(i)        = projected coord 
c
c  Note:  rlond(i) & u(i) can safely be same array, as can be rlatd(i) & v(i).
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlon0
c     real sgn
c     real factps
c     integer n
      dimension rlatd(n),rlond(n)
      dimension u(n),v(n)
c
c
c  Define degrees to radians factor symbolic constant
c
      parameter(DEG2RAD=.0174532925199433)
      parameter(POLE=89.99)
c
c
c  Compute projected coordinates for each input longitude & latitude
c
      do 3100 i=1,n
       phi = max( -POLE, min( POLE, rlatd(i) ) )
       r = factps * tan( DEG2RAD * .5 * ( 90. + sgn * phi ) )
       theta = DEG2RAD * sgn * ( 90. - ( rlond(i) - rlon0 ) )
       u(i) = r * cos(theta)
       v(i) = r * sin(theta)
 3100 continue
c
c
c  Return to caller with u(i),v(i)
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine invplc(cone,rlon0,sgn,n, u,v, rlond,rlatd)
c
c
c  @(#) invplc.f  McKie Jul-1989
c  This routine computes the inverse transform from projection space back to
c  lat,lon space, given a particular point in projection space, for the
c  lambert conformal projection with 2 parallels.
c
c  Input:
c                cone = cone constant
c               rlon0 = central longitude
c                 sgn = -1. for northern hemisphere, +1. for southern
c                   n = number of coordinates to be inverse projected
c           u(i),v(i) = projected coordinates
c
c  Output:
c     rlatd(i),rlond(i) = latitude, longitude coord corresponding to u(i),v(i)
c
c  Note:  rlond,rlatd can safely be the same array as u,v in calling routine.
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real cone
c     real rlon0
c     real sgn
c     integer n
      dimension u(n),v(n)
      dimension rlatd(n),rlond(n)
c
c
c  Define local constants
c
      parameter(RAD2DEG=57.2957795130823)
c
c
c  Scan the list of projected space u,v; & transform each pair to lon,lat coord
c
      factor = 2. * sgn
      expon = 1. / cone
      do 3100 i=1,n
       r = sqrt( u(i)**2 + v(i)**2 )
       rlond(i) = rlon0 + RAD2DEG * atan2(u(i),sgn*v(i)) / cone
       rlatd(i) = factor * ( atan( r**expon ) * RAD2DEG - 45. )
 3100 continue
c
c
c  Return to caller with inverse projection coord in rlatd(i),rlond(i)
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine invpme(rlon0,factme, n, u,v, rlond,rlatd)
c
c
c  @(#) invpme.f  McKie Mar-1992
c  This routine computes the inverse transform from projection space back to
c  lon,lat space, given a particular point (u,v) in projection space, for the
c  Mercator projection.  [See Hzliner & Williams: Num Pred & Dyn Met, p 13.]
c
c  Input:
c               rlon0 = central longitude
c              factme = Mercator factor, cos(rlat0)
c                   n = number of coordinates to transform
c           u(i),v(i) = projected coordinates
c
c  Output:
c     rlatd(i),rlond(i) = latitude, longitude coord corresponding to u(i),v(i)
c
c  Note:  rlond,rlatd can safely be the same array as u,v in calling routine.
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlon0
c     real factme
c     integer n
      dimension u(n),v(n)
      dimension rlatd(n),rlond(n)
c
c
c  Define local constants
c
      parameter(RAD2DEG=57.2957795130823)
      parameter(ANGLEREF=-180.)
c
c
c  Scan the list of projected space u,v; & transform each pair to lon,lat coord
c
      do i=1,n

       angle = rlon0 + RAD2DEG * ( u(i) / factme )
       call angrng( angle, ANGLEREF, rlond(i) )

c      theta = atan( exp( v(i) / factme ) )                     ! slower method
c      rlatd(i) = RAD2DEG * asin( 2. * ( sin(theta)**2 ) - 1. ) ! slower method

       r2 = exp( v(i) / factme ) ** 2                           ! faster method
       rlatd(i) = RAD2DEG * asin( ( r2 - 1 ) / ( r2 + 1 ) )     ! faster method

      enddo
c
c
c  Return to caller with inverse projection coord in rlatd(i),rlond(i)
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      subroutine invpps(rlon0,sgn,factps, n, u,v, rlond,rlatd)
c
c
c  @(#) invpps.f  McKie Sep-1991
c  This routine computes the inverse transform from projection space back to
c  lon,lat space, given a particular point (u,v) in projection space, for the
c  polar stereographic projection.
c
c  Input:
c               rlon0 = central longitude
c                 sgn = -1. for northern hemisphere, +1. for southern
c              factps = polar stereographic factor, abs(sin(phi0) - sgn)
c                   n = number of coordinates to transform
c           u(i),v(i) = projected coordinates
c
c  Output:
c     rlatd(i),rlond(i) = latitude, longitude coord corresponding to u(i),v(i)
c
c  Note:  rlond,rlatd can safely be the same array as u,v in calling routine.
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
c  Declare subprogram args
c
c     real rlon0
c     real sgn
c     real factps
c     integer n
      dimension u(n),v(n)
      dimension rlatd(n),rlond(n)
c
c
c  Define local constants
c
      parameter(RAD2DEG=57.2957795130823)
      parameter(ANGLEREF=-180.)
c
c
c  Scan the list of projected space u,v; & transform each pair to lon,lat coord
c
      do 3100 i=1,n
       r = sqrt( u(i)**2 + v(i)**2 )
       angle = 90. + rlon0 - RAD2DEG * sgn * atan2( v(i), u(i) )
       call angrng( angle, ANGLEREF, rlond(i) )
       rlatd(i) = ( 2. * RAD2DEG * atan( r / factps ) - 90. ) * sgn
 3100 continue
c
c
c  Return to caller with inverse projection coord in rlatd(i),rlond(i)
c
      return
      end
c
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c
c
      SUBROUTINE ANGRNG(A,AREF,ANEW)
C
C
C  THIS ROUTINE COMPUTES AN EQUIVALENT ANGLE TO A SPECIFIED ANGLE
C  WHICH LIES WITHIN A SPECIFIED 360 DEGREE RANGE OF VALUES.
C
C  INPUT:      A = ORIGINAL ANGLE (DEGREES)
C           AREF = REFERENCE ANGLE FOR RANGE [AREF,AREF+360.]
C
C  OUTPUT:  ANEW = ANGLE EQUIVALENT TO A WHICH LIES WITHIN [AREF,AREF+360]
C
C  IT IS DESIGNED TO RUN ON THE CRI PDP-11/70, RSX11M, F4P.
C
C  @(#) angrng.f  McKie  Jun-1981
c
c
c  Use implicit variable typing header file
c
      include 'precision.h'
c
c
      R = (A - AREF) / 360.
      IF(R .LT. 0.) R = R - 1.
      ANEW = A - int(R) * 360.
      RETURN
C
      END
