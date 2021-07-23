subroutine readinput()
   use comvar
   implicit none
  
   open(unit=20, file="param.in", status="old", action="read" )
   read(20,*) gfilename  !File name of grid
   gfilename = trim(gfilename)
   !gfilename = "feflo.domn.naca0012"
   !gfilename = trim(gfilename)
   write(*,*) "The grid file is ", gfilename
   read(20,*) 
   read(20,*) case1
   write(*,*) case1
   read(20,*) cfl
   write(*,*) cfl

   read(20,*) itmax
   write(*,*) itmax
   read(20,*) itsave
   write(*,*) itsave
   read(20,*) restol
   write(*,*) restol
   read(20,*) limtype
   write(*,*) limtype
   read(20,*) fluxtype
   write(*,*) fluxtype
   read(20,*) lim_par
   write(*,*) lim_par

   close(20)

end subroutine readinput
