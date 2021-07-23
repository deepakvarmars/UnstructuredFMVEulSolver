!Initial condition for cylinder test case
subroutine userinit()
   use comvar
   implicit none
   real    :: aoa, m, ss, velmag

   write(*,*) "Initial Values of the Euler variables"

   aoa = 0.0  !Angle of attack in degrees

   aoa = aoa*c_pi/180.0

   m = 0.38 !Mach no
   velmag = 1.0 !Magnitude of velocity

   ss = velmag/m

   velin(1) = velmag*cos(aoa) 
   velin(2) = velmag*sin(aoa) 

   rhoin = 1.0

   gam = 1.4

   prein = rhoin*ss**2/gam

end subroutine userinit

subroutine useric(q)
   use comvar
   implicit none
   real, intent(inout) :: q(4)

   q(1) = rhoin
   q(2) = velin(1)
   q(3) = velin(2)
   q(4) = prein
   
end subroutine useric
