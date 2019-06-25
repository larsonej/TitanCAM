;
;...Numerical solution for coagulation given by Golovin [1963]
;

f = ' '
x = (findgen(40001)+1.)/1000.
phi = exp(-x)

plot_io, x, phi, yrange=[1.e-10,1.]
print,total(phi)
read,f

for itime = 1, 10 do begin

  t = itime*0.01
  tau = 1. - exp(-t)
  phi = (1.-tau) * exp(-x*(tau+1.)) / (x*sqrt(tau)) * beseli(2.*x*sqrt(tau), 1)

  oplot, x, phi
  print,total(phi)
read,f

endfor

stop
end
