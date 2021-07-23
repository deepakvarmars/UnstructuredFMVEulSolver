!Subroutine to set the initial conditions of primitive variables
subroutine setic()
   use comvar
   implicit none
   integer :: ieqns, ielem

   call userinit()

   unkel = 0.0
   cunkel = 0.0
   write(*,*) "Set Initial Condition"
   do ielem =1 ,nelem
      call useric(unkel(:,ielem))
      !call primtocon(unkel(:,ielem), cunkel(:,ielem))
      cunkel(1,ielem) = unkel(1,ielem)
      cunkel(2,ielem) = unkel(2,ielem)*unkel(1,ielem)
      cunkel(3,ielem) = unkel(3,ielem)*unkel(1,ielem)
      cunkel(4,ielem) = unkel(4,ielem)/(gam - 1.0) + &
                        0.5*unkel(1,ielem)*(unkel(2,ielem)**2 + &
                        unkel(3,ielem)**2)
      do ieqns = 1, neqns
         if(isnan(unkel(ieqns,ielem))) then
            write(*,*) "Set IC Nan"
            write(*,*) "eqn no, elem no ",ieqns,ielem
            call exit(0)
         endif
      enddo
   enddo

 
end subroutine setic

