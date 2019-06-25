;
;...Plot results from coag simulation
;

title='aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
simtitle=title

close,1
openr,1,'model_his.out', /f77_unformat

readu,1, title
readu,1, simtitle
;print,title
;print,simtitle

nx = 1L  &  ny = 1L  &  nz = 1L  &  nbins = 1L  &  nelem = 1L  &  ngroups = 1L

readu,1, nx, ny, nz, nbins, nelem, ngroups

;print,nx,ny,nz

r = dblarr(ngroups,nbins)
dr = r  &  rmass = r  &  dm = r
c = dblarr(nbins,nelem,50)
totvol = dblarr(50)
cr = dblarr(nbins,nelem)

readu,1, r, dr, rmass, dm

itime = 1L
time = 0.d0

for i = 0,49 do begin

   readu,1, itime, time
;   print,itime, time

   readu,1, cr
   c(0,0,i) = cr(0:nbins-1,0:nelem-1)
   totvol(i) = total(cr(*,1)*r(1,*)^3 + cr(*,0)*r(0,*)^3)
;   print, i, totvol(i)

endfor

;stop

plot_oo, r(1,*), c(*,1,0), yrange=[1.e-3,1.e3], psym=-1
for i = 0,49 do oplot, r(0,*), c(*,0,i), psym=-1
for i = 0,49 do oplot, r(1,*), c(*,1,i), psym=-1, linestyle=2


stop
end
