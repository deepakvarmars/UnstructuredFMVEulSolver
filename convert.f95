!Convert primitive to conservative variables
subroutine primtocon(q,u)
   use comvar
   implicit none
   real(8), intent(in)   :: q(neqns)
   real(8), intent(out)  :: u(neqns)

   u(1) = q(1) 
   u(2) = q(2)*q(1) 
   u(3) = q(3)*q(1) 
   print*, q(1), q(3), u(3)
   u(4) = q(4)/(gam - 1.0d0) + 0.5d0*q(1)*(q(2)**2 + q(3)**2)

   print*, q
   print*, u
   stop
   return

end subroutine primtocon

!Convert conservative to primitive variables
subroutine contoprim(u,q)
   use comvar
   implicit none
   real(8), intent(in)   :: u(neqns)
   real(8), intent(out)  :: q(neqns)

   q(1) = u(1) 
   q(2) = u(2)/u(1) 
   q(3) = u(3)/u(1) 
   q(4) = (gam - 1.0d0)*(u(4) - 0.5d0*q(1)*(q(2)**2 + q(3)**2))
   q(4) = max(1.0d-7,q(4))

   return

end subroutine contoprim
