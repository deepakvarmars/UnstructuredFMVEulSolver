subroutine preprocess()
   use comvar
   implicit none
   
   write(*,*) "Pre-processing"
   !Find elements surrounding points
   write(*,*) "Find the elements surrounding points"
   call elemsupt()


   !Find elements surrouding elements
   write(*,*) "Find the elements surrounding elements"
   call elemsuelem()

   !Find the elements containing bfaces and no.of interior faces
   write(*,*) "Find the elements containing boundaries and no. of internal faces"
   call nfaces()

end subroutine preprocess

!Element surrounding points
subroutine elemsupt()
   use comvar
   implicit none
   integer :: inode, ielem, ipoin, istor
   
   esup2 = 0

   !Element pass - 1
   do ielem = 1, nelem
      do inode = 1, nnode
         esup2(inpoel(inode,ielem) + 1) = esup2(inpoel(inode, ielem)+1) + 1
      enddo
   enddo
   
   !Storage/Reshuffling pass - 1
   do ipoin = 2, npoin + 1
      esup2(ipoin) = esup2(ipoin) + esup2(ipoin - 1)
   enddo

   nesup = esup2(npoin+1)
   allocate(esup1(nesup)) 
   esup1 = 0

   !Element pass - 2
   do ielem = 1, nelem
      do inode = 1, nnode
         ipoin = inpoel(inode,ielem)
         istor = esup2(ipoin) + 1
         esup2(ipoin) = istor
         esup1(istor) = ielem
      enddo
   enddo

   !Storage/Reshuffling pass - 2
   do ipoin = npoin + 1, 2, -1
      esup2(ipoin) = esup2(ipoin-1)
   enddo

   esup2(1) = 0

end subroutine elemsupt


!Elements surrounding elements
subroutine elemsuelem()
   use comvar
   implicit none
   integer :: ielem, ifael
   integer :: ip1, ip2
   integer :: istor, jelem, icoun, ieadj
   integer,allocatable,dimension(:) :: lpoin
   integer, dimension(2,3) :: lhelp

   allocate(lpoin(npoin))

   lhelp = reshape((/2,3,  3,1,  1,2/),shape(lhelp))

   lpoin = 0
   esuel = 0

    do ielem = 1,nelem
      do ifael = 1,nfael

        ieadj = 0

        ! obtain nodes of this face
        ip1 = inpoel(lhelp(1,ifael),ielem)
        ip2 = inpoel(lhelp(2,ifael),ielem)
        lpoin(ip1) = 1    ! mark in lpoin
        lpoin(ip2) = 1

        ! loop over the elements surrounding point ip1
        do istor = esup2(ip1)+1,esup2(ip1+1)

          jelem = esup1(istor)  ! element number

          if(jelem .ne. ielem) then

            icoun = lpoin(inpoel(1,jelem)) + lpoin(inpoel(2,jelem)) &
              + lpoin(inpoel(nnode,jelem))

            if(icoun == nnode-1) then

              ieadj = jelem
              go to 10

            endif

          endif

        enddo

        10 continue  ! the adjacent element is found

        esuel(ifael,ielem) = ieadj
        lpoin(ip1) = 0
        lpoin(ip2) = 0

      enddo  ! over the faces
    enddo  ! over the elements

end subroutine elemsuelem




!Find number of faces and boundary elements : nface 
subroutine nfaces()
   use comvar
   implicit none
   integer :: ifac, ielem, inode, ip1, ip2
   integer :: fcn1, fcn2, jelem
   integer, allocatable, dimension(:,:) :: intfactmp

   allocate(intfactmp(4,mafac))
   !Find the elements contaning boundary faces
   
   !Loop over the boundary faces
   do ifac = 1, nbfac

      !Store the bface nodes
      ip1 = bface(1,ifac)
      ip2 = bface(2,ifac)

      do ielem = 1, nelem    ! loop over elems for finding b_elem
         do inode = 1, nnode ! loop over the vertices of the elems
          
           fcn1  = mod(inode,nnode) + 1   ! first node of face
           fcn2  = mod(fcn1,nnode) + 1 ! second node of face

           if(inpoel(fcn1, ielem) == ip1 .and. inpoel(fcn2, ielem) == ip2) then  !check if nodes correspond to bface
              intfactmp(1,ifac) = ielem          ! Boundary elements
              intfactmp(2,ifac) = nelem + ifac   ! Ghost element 
              intfactmp(3,ifac) = ip1            ! First node of bface
              intfactmp(4,ifac) = ip2            ! Second node of bface
           endif

         enddo
      enddo
   enddo

   nafac = 0 !Initialise the number of face(inter + bound)
   nafac = nbfac 
   !Find the number of internal faces
   do ielem = 1, nelem      !Loop over the elements
      do inode = 1, nnode   !Loop over each node of the element
  
         fcn1 =  mod(inode,nnode) + 1
         fcn2 =  mod(fcn1,nnode) + 1

         jelem = esuel(inode,ielem)

         if(jelem > ielem .and. jelem <= nelem) then    !Make sure its not a boundary face

            nafac = nafac + 1

            if(nafac > mafac+1) then
               write(*,*) "Err : Increase mafac"
               call exit(0)
            endif

            intfactmp(1,nafac) = ielem
            intfactmp(2,nafac) = jelem
            intfactmp(3,nafac) = inpoel(fcn1,ielem)
            intfactmp(4,nafac) = inpoel(fcn2,ielem)

         endif

      enddo

   enddo 

   write(*,*) "Number of all faces = ", nafac

   allocate(intfac(4,nafac))
 
   do ifac=1,nafac
      intfac(1:4, ifac) = intfactmp(1:4,ifac)
   enddo 
  
   deallocate(intfactmp)

   return

end subroutine nfaces
