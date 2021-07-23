subroutine timestep()
   use comvar
   implicit none
   integer :: inode, ielem
   integer :: inode1, inode2, ip1, ip2, jelem
   real :: p1, p2, r1, r2, c1, c2, area, vel1, vel2
   real :: nx, ny

   wgt = 0.0
   dt1 = 0.0
   do ielem = 1, nelem
      
      do inode=1,3
      
         inode1 = mod(inode,nnode) + 1
         inode2 = mod(inode1,nnode) + 1
       
         jelem = esuel(inode,ielem)
        
         !Bugs for cylinder mesh --- jelem = 0
         if(jelem == 0) cycle
            
         ip1 = inpoel(inode1,ielem)
         ip2 = inpoel(inode2,ielem)

         nx = coord(2,ip2) - coord(2,ip1) 
         ny = coord(1,ip1) - coord(1,ip2)
         area = sqrt(nx**2 + ny**2)
         nx = nx/area 
         ny = ny/area 

         r1 = max(1d-6,unkel(1,ielem))
         p1 = max(1d-6,unkel(4,ielem))
         c1 = sqrt(gam*p1/r1)
         r2 = max(1d-6,unkel(1,jelem))
         p2 = max(1d-6,unkel(4,jelem))
         c2 = sqrt(gam*p2/r2)

         vel1 = unkel(2,ielem) + unkel(2,jelem)
         vel2 = unkel(3,ielem) + unkel(3,jelem)

         wgt(ielem) = wgt(ielem) + 0.5d0*area*(c1 + c2 + abs(vel1*nx + vel2*ny))
      enddo
      dt1(ielem) = cfl*geoel(3,ielem)/wgt(ielem)
   enddo

   dt = minval(dt1)


end subroutine timestep
