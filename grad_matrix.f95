!Generate the left hand side matrix of least squares approach for computing gradients

subroutine grad_matrix_generate()
   use comvar
   implicit none
   integer :: ielem, jelem, ifael, ifac
   real    :: dxj, dyj, wj

   wmat = 0.0d0

   do ielem = 1, nelem   !Loop over each element

      do ifael= 1, nfael  !Loop over faces of the element

         jelem = esuel(ifael, ielem)

         if(jelem > nelem) then !This is a boundary face
            ifac = jelem - nelem
            
            dxj = geofa(1,ifac)
            dyj = geofa(2,ifac)

            wj = 1.0d0/geofa(3,ifac)**2

         else
 
            dxj = geoel(1,jelem) - geoel(1,ielem)
            dyj = geoel(2,jelem) - geoel(2,ielem)

            wj =  1.0d0/(dxj**2 + dyj**2)
         endif

         wmat(1,ielem) = wmat(1,ielem) + wj*dxj**2
         wmat(2,ielem) = wmat(2,ielem) + wj*dxj*dyj
         wmat(3,ielem) = wmat(3,ielem) + wj*dyj**2

         
      enddo

   enddo
   
   return

end subroutine grad_matrix_generate
