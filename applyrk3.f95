subroutine applyrk3(rk)
   use comvar
   implicit none
   integer, intent(in) :: rk
   integer :: ielem
   real(8) :: brk

   brk = 1.0d0 - ark(rk)

   do ielem = 1, nelem
      cunkel(:,ielem) = ark(rk)*q0(:,ielem) + brk*(cunkel(:,ielem) &
                        + dt*rhs(:,ielem)/geoel(3,ielem))
      unkel(1,ielem) = cunkel(1,ielem)
      unkel(2,ielem) = cunkel(2,ielem)/cunkel(1,ielem)
      unkel(3,ielem) = cunkel(3,ielem)/cunkel(1,ielem)
      unkel(4,ielem) = (gam - 1.0)*(cunkel(4,ielem) - &
                        0.5*unkel(1,ielem)*(unkel(2,ielem)**2 +&
                         unkel(3,ielem)**2))
   enddo

   return

end subroutine applyrk3
