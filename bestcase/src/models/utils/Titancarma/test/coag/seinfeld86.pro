;
;...Plot coagulation model results for self-preserving distribution [Seinfeld, 1986]
;

;
;
;...Read history file from AERMOD coagulation simulation
;
;
title='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
simtitle=title
 
close,1
openr,1,'model_his.out', /f77_unformat
 
;...Read titles

readu,1, title
readu,1, simtitle
print,title
print,simtitle
 
;...Read array dimensions
 
nx = 1L  &  ny = 1L  &  nz = 1L  &  nbins = 1L  &  nelem = 1L  &  ngroups = 1L
 
readu,1, nx, ny, nz, nbins, nelem, ngroups
 
print,nx,ny,nz
 
;...Read radius, mass, and particle number concentrations
 
r = dblarr(nbins,ngroups)
dr = r  &  rmass = r  &  dm = r
c = dblarr(nbins,nelem,101)
totvol = dblarr(101)
cr = dblarr(nbins,nelem)
 
readu,1, r, dr, rmass, dm
 
itime = 1L
time = 0.d0
 
for i = 0,100 do begin
 
   readu,1, itime, time
 
   readu,1, cr
   c(0,0,i) = cr(0:nbins-1,0:nelem-1)
   totvol(i) = total(cr(*,1)*r(1,*)^3 + cr(*,0)*r(0,*)^3)
;   print, i, totvol(i)
 
endfor
 
 
tau = findgen(101)/10.

;...Analytic solution

m = fltarr(40)
m(0) = 1.e-21
mrat = 2.
for i = 1,39 do m(i) = m(i-1)*mrat
delm = 2*m * (mrat-1.) / (mrat+1.)
tau = findgen(101)/10.

n0 = 1.
i0 = 9

nana = fltarr(40,101)

for itau = 0,100 do nana(0,itau) = n0 * reform(dm(0:nbins-1,0)) / rmass(i0,0) / (1.+tau(itau))^2 *     $
                                     exp( -reform(rmass(0:nbins-1,0))/rmass(i0,0)/(1.+tau(itau)) )

;...Plot size distributions

psinit,0
device,filename='seinfeld86.ps'

!p.position = [0.15,0.07,0.9,0.7]
plot_oo, rmass(*,0), c(*,0,0)*4.2e-21/dm(*,0), yrange=[1.e-6,0.4], xrange=[4.2e-21,1.e-14], psym=-1, /noerase, /nodat, $
         xtitle='Mass', ytitle='dN/dlog(m)', title='Seinfeld [1985] self-preserving solution'
for i = 00,20,20 do oplot, rmass(*,0), c(*,0,i)*rmass(*,0)/dm(*,0), psym=-1
for i = 100,100,100 do oplot, rmass(*,0), c(*,0,i)*rmass(*,0)/dm(*,0), psym=-1

for i = 00,20,20 do oplot, rmass(*,0), nana(*,i)*rmass(*,0)/dm(*,0), psym=-1, linestyle=2
for i = 100,100,100 do oplot, rmass(*,0), nana(*,i)*rmass(*,0)/dm(*,0), psym=-1, linestyle=2

xline = [3.e-17,1.e-16]
yline = [0.2,0.2]  &  oplot,xline,yline  &  xyouts, xline(1)+0.1e-16,yline(0)*0.93, 'Numerical: Mrat='+ $
                                                    string(rmass(1,0)/rmass(0,0),form='(F5.2)'), size=1.05
yline = [0.09,0.09]  &  oplot,xline,yline,linestyle=2  &  xyouts, xline(1)+0.1e-16,yline(0)*0.93, 'Analytic', size=1.05

device,/close
spawn,"gs seinfeld86.ps"

;plot_oo, m, n, xrange=[1.e-21,1.e-11], yrange=[1.e-8,1.]

stop
end
