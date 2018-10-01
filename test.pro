function fkernel,u

   fac = 10.666666666667 + u * u * (32.0 * u - 38.4)
   jj=where(u GE 0.5) 
  IF jj(0) GT -1 THEN  fac(jj) = 21.333333333333 - 48.0 * u(jj) + $
                              38.4 * u(jj) * u(jj) - 10.666666666667 * u(jj) * u(jj) * u(jj) - $
                              0.066666666667 / (u(jj) * u(jj) * u(jj))

   return,fac
end

function fkerneli,u

   fac = 10.666666666667 + u * u * (32.0 * u - 38.4)
   jj=where(u GE 0.5) 
  IF jj(0) GT -1 THEN  fac(jj) = 21.333333333333 - 48.0 * u(jj) + $
                              38.4 * u(jj) * u(jj) - 10.666666666667 * u(jj) * u(jj) * u(jj)

   return,fac
end

pro test
  N=100
  u=findgen(N)/(N-1)
  k=fkernel(u)
  ki=fkerneli(u)

  plot,u,k
  oplot,u,ki,col=mycolor(4),l=2

  stop
end