! 2nd Order FV Solver for Euler equations
! Author - Deepak Varma Ravivarma Sylaja
! NC State ID - 200281115
program main
   use comvar
   implicit none
    
   write(*,*) "Read input file" 
   !Read input file
   call readinput()

   write(*,*) "Read grid file" 
   !Read the grid file
   call readgrid()

   !Pre-processing
   call preprocess()


   ! Find geometric values of mesh (i.e. area of elems, normals etc)
   call geombuild()

   !LHS of least square approximation build gradients
   if(limtype > 1) then
      call grad_matrix_generate()
   endif

   ark(1) = 0.0d0
   ark(2) = 0.75d0
   ark(3) = 1.0d0/3.0d0

   !Processing
   call runsim()

   write(*,*) "Finished simulation"

end program main
