subroutine readgrid()
   use comvar
   implicit none
   integer   ::  i,junkint
   real      ::  junkrl
 
   write(*,*) "Read the grid file"
   open(unit = 10, file = gfilename, status = "OLD", action="READ")
   read(10,*) junkint
   do i = 1, junkint+1
      read(10,*) 
   end do
   read(10,*) ndimn, nnode
   read(10,*) 
   read(10,*) nelem, npoin, nbfac, junkrl

   write(*,*) "No.of elements = ", nelem
   write(*,*) "No.of points = ", npoin
   write(*,*) "No.of boundary faces = ", nbfac
   nconi = ndimn + 2 
   write(*,*) "Allocating the necessary arrays"
   call allocator()
   write(*,*) "Finished allocating"

   write(*,*) "Read the connectivity matrix"
   !Read in connectivity matrix
   read(10,*)
   do i=1,nelem
      read(10,*) junkint, inpoel(1,i), inpoel(2,i), inpoel(3,i)
   end do
   read(10,*)
  !read the co-ordinates of each global node
   write(*,*) "Read the co-ordinates of each node"
   do i = 1, npoin
      read(10,*) junkint, coord(1,i), coord(2,i)
   end do
   do i = 1, npoin+1
      read(10,*) 
   enddo 
   !Read data for the boundary faces
   read(10,*)
   write(*,*) "Read boundary face data"
   do i = 1, nbfac
      read(10,*)  junkint, bface(1,i), bface(2,i), bface(3,i), junkint, bface(4,i)
   end do
   write(*,*) "Finished reading grid file" 
   close(10)

end subroutine  readgrid
