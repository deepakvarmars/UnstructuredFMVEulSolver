subroutine runsim()

   use comvar
   implicit none
   integer :: j, k, it, rk
   integer :: solfile
   real    :: resnorm0, resnorm
   
   !Set the initial condition
   call setic()

   solfile = 0
   resnorm = 1.0e20
   it = 0
   time = 0.0
   !Write the initial solution
   call writesol(solfile)

   do while(resnorm > restol .and. it < itmax)

      !Find the timestep for the iteration
      call timestep()
      !HARDCODING SMALL TIME STEP
      !dt = 10e-7

      !Store older conservative variables
      q0 = cunkel
      
      !Use RK3 scheme
      do rk = 1,3

         !Set boundary condition
         call setbc()

         !Find the numerical fluxes and compute RHS
         call computerhs()

         !Use time-integration technique
         call applyrk3(rk)

      enddo

      !Find the l2norm of relative value of rhs/residual
      resnorm = 0.0
      do k = 1, nelem
         do j = 1, neqns
            resnorm = resnorm + rhs(j,k)**2
         enddo
      enddo
      resnorm = sqrt(resnorm/(neqns*nelem))
      if(isnan(resnorm)) then
         write(*,*) "Err: RHS NaN"
         write(*,*) "iteration no = ", it
         call exit(0)
      endif

      if(it == 0) then
         resnorm0 = resnorm
         write(*,*) "Initial residual norm =", resnorm0
         write(*,*) "it   dt	   time	   resnorm"
      endif
      if(resnorm0 > 10e-13) resnorm = resnorm/resnorm0 !Compute relative residual
      
      time = time + dt
      it = it + 1
 
      !Output the relative residuals
      write(*,'(i8,3e14.5)') it, dt, time, resnorm

      !Save solutions
       if(mod(it,itsave)==0 .or. resnorm < restol) then
          solfile = solfile + 1
          call writesol(solfile)
       endif

   enddo
   
   !deallocate(q0)

end subroutine runsim
