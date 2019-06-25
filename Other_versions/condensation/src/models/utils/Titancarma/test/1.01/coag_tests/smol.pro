;
;...Plot results from coag simulation
;
 
title='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
simtitle=title
 
close,1
openr,1,'model_his.out', /f77_unformat
 
readu,1, title
readu,1, simtitle
print,title
print,simtitle
 
nx = 1L  &  ny = 1L  &  nz = 1L  &  nbins = 1L  &  nelem = 1L  &  ngroups = 1L
 
readu,1, nx, ny, nz, nbins, nelem, ngroups
 
print,nx,ny,nz
 
r = dblarr(ngroups,nbins)
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
totc = fltarr(101)
for i = 0,100 do totc(i) = total(c(*,0,i))

psinit,0
device,filename='smol.ps'

;...Plot total number concentration and number in first bin

!p.position = [0.15,0.57,0.9,1.]
plot_io, tau, c(0,0,*), xrange=[0,10], yrange=[0.005,1.], psym=0, xtitle='tau', ytitle='Number Concentration',  $
                        tit='Smoluchowski solution'
oplot, tau, totc, psym=0

;...Analytic solution

c1 = 1. / (1.+tau)^2
ctot = 1. / (1.+tau)
oplot, tau, c1, linestyle=2
oplot, tau, ctot, linestyle=2

xline = [5.,6.]
yline = [0.7,0.7]  &  oplot,xline,yline  &  xyouts, xline(1)+0.1,yline(0)*0.93, 'Numerical: Mrat='+ $
                                                    string(rmass(0,1)/rmass(0,0),form='(F5.2)'), size=1.1
yline = [0.5,0.5]  &  oplot,xline,yline,linestyle=2  &  xyouts, xline(1)+0.1,yline(0)*0.93, 'Analytic', size=1.1
xyouts,5.,0.3, 'Top curves: Total Number', size=1.1
xyouts,5.,0.22, 'Bottom curves: Number in bin 1', size=1.1

;...Plot size distributions

!p.position = [0.15,0.07,0.9,0.5]
plot_oo, rmass(0,*), c(*,0,0)*4.2e-21/dm(0,*), yrange=[1.e-3,0.4], xrange=[4.2e-21,1.e-19], psym=-1, /noerase, /nodat, $
         xtitle='Mass', ytitle='dN/dm'
for i = 10,20,10 do oplot, rmass(0,*), c(*,0,i)*4.2e-21/dm(0,*), psym=-1
for i = 50,100,50 do oplot, rmass(0,*), c(*,0,i)*4.2e-21/dm(0,*), psym=-1

;...Analytic solution

mana = (1.+findgen(50))*4.2e-21
k50 = findgen(50)
nana = fltarr(50,101)
for itau = 0,100 do nana(0,itau) = tau(itau)^(k50(0:49)) / (1.+tau(itau))^(k50(0:49)+2)

for i = 10,20,10 do oplot, mana, nana(*,i), psym=-1, linestyle=2
for i = 50,100,50 do oplot, mana, nana(*,i), psym=-1, linestyle=2

device,/close
spawn,"gs smol.ps"

stop
end
