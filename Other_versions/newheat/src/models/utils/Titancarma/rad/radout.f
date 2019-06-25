      subroutine radout
c
C     **************************************************************
C     *  Purpose             :  Prints out radiative transfer      *
C     *                         results.                           *
C     **************************************************************

      include 'globrad.h'

      write(LUNORAD,*) 'radout: u0 = ',u0
      write(LUNORAD,*) 'radout: alb(toa)=',alb_toai
      write(LUNORAD,*) 'radout: alb(tom)=',alb_tomi
      write(LUNORAD,*) 'solar energy absorbed=',fnetbs(nlayer)-fnetbs(1)
      write(LUNORAD,*) 'ir energy absorbed=',fnetbi(nlayer)-fnetbi(1)
C
C       PRINT OUT THE INPUT VARIABLES
C
        WRITE(LUNOPRT,560)
 560    FORMAT(" RADOUT:",/,
     1         " j   p(j)    press(j)    t(j) ",
     2         "     tt(j)  rdh2o(j)  ctot       firu     fird  ",
     3         "    fsLu     fsLd   ")
C
        DO 565 J = 1, NVERT
           ctot = 0.
           do ig = 1, NGROUP
             do i = 1, NRAD
               ctot = ctot + caer(i,j,ig)
             enddo
           enddo
           WRITE(LUNOPRT,562) J,P(J),PRESS(J),T(J),TT(J),
     1                  RDH2O(J),Ctot,fupbi(j),fdownbi(j),
     2                  fupbs(j),fdownbs(j)
 562       FORMAT(I3,11(1PE9.2))
 565    CONTINUE
        WRITE(LUNOPRT,563) NLAYER,PRESS(NLAYER),TT(NLAYER)
 563    FORMAT(I3,9X,1PE9.2,9X,1PE9.2)
C
C       PRINT OUT THE CALCULATED VARIABLES
C
        WRITE(LUNOPRT,526)
 526    FORMAT("***delta scaled***:",/
     1     "  solnet     xirdown     irup       u0",
     2     "       opd(1)     opd(8)     opd(16)    opd(26)    opd(38)")
        WRITE(LUNOPRT,530) SOLNET,XIRDOWN,XIRUP,
     1         U0,OPD(1,NLAYER),OPD(8,NLAYER),OPD(16,NLAYER),
     2         OPD(26,NLAYER),OPD(38,NLAYER)
 530    FORMAT(9(1PE10.3,1X))
C
        WRITE(LUNOPRT,536)
 536    FORMAT("***delta scaled***:"/,
     1         "              j    heats",
     1    "       heati        heat         wol         gol        opd",
     2    "        taul")
C
        DO 550 J         =  1,NVERT
           WRITE(LUNOPRT,540) J,HEATS(J),HEATI(J),
     1                  HEAT(J),WOL(8,J),GOL(8,J),OPD(8,J),TAUL(8,J)
 540       FORMAT(10X,I4,7(1PE11.3,1X))
 550    CONTINUE

      i_more_outpt = 0
      if (i_more_outpt.ne.0) then
        j = NVERT/2. + 1.5
        write(LUNORAD,780) j-1
 780    format(//,'RADOUT: unscaled optical properties for layer',i3,/,
     1  '                ----------- cloud --------------- ',
     2  ' ----------cloud and gas---------',/
     3  '   i     wave        tau         w0         g0 ',
     4  '       tau         w0         g0       ')
        do 790 i = 1, nwave
          index = nprobi(i,1)
          w0cloud = taus(i,j)/taua(i,j)
          write(LUNORAD,781) i,wave(i),taua(i,j),w0cloud,goL(index,j),
     1               utauL(index,j),uw0(index,j),ug0(index,j)
 781      format(1x,i4,7(1pe11.2))
          do 790 ii = 2, nprobi(i,2)
            index = index + 1
            write(LUNORAD,784) utauL(index,j),uw0(index,j),ug0(index,j)
 784        format(49x,3(1pe11.2))
 790    continue
        write(LUNORAD,785) nwave+1,wave(nwave+1)
 785    format(1x,i4,1pe11.2)
      endif

      write(LUNORAD,566)
 566  format(//,' RADOUT: Top of atmosphere radiative fluxes:',/
     1   '       -------------------solar----------------- ',
     2   ' --infra-red (top of model)--'/,
     3   '    i    wave         up       down     albedo   ',
     4   '       wave     up')

      write(LUNORAD,567) (i,wave(i),fsLu(i),fsLd(i),alb_toa(i),
     1  .5*(wave(i+nsol)+wave(i+1+nsol)),firu(i),i=1,nir)
 567  format((1x,i4,3(1pe11.2),0pf9.3,2x,2(1pe11.2)))
      write(LUNORAD,568) (i,wave(i),fsLu(i),fsLd(i),alb_toa(i),
     1            i=nir+1,nsol)
 568  format ((1x,i4,3(1pe11.2),0pf9.3))

      write(LUNORAD,782) tslu,tsld,alb_toai,tiru
 782  format(' totals:        ',2(1pe11.2),0pf9.3,11x,1pe11.2)

      write(LUNORAD,556)
 556  format(//,' RADOUT: total unscaled optical depth',/
     1   '    i    wave        opd')
      write(LUNORAD,537) (i,wave(i),uopd(i,nlayer),i=1,nwave)
 537  format((1x,i4,1p2e10.2,0p))

      return
      end
