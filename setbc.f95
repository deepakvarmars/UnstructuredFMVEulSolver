subroutine setbc()
   use comvar
   implicit none
   integer :: ifac, ie, ge
   real    :: nx, ny, lng !unit normals and lng
   real    :: c0, s0, c1, s1
   real    :: vn0, vt0, vt1
   real    :: cg, vng, vtg 
   real    :: vn1, mn1

   !Here 1  as postfix : internal   ; 0 as postfix : farfield/initial values

   do ifac = 1, nbfac
  
      lng = geofa(3,ifac)  
      nx  = geofa(1,ifac)/lng
      ny  = geofa(2,ifac)/lng

      ie  = intfac(1,ifac)
      ge  = intfac(2,ifac)

      if(bface(3,ifac) == 2) then    !Slip-wall
   
         unkel(1,ge) = unkel(1,ie)
         unkel(2,ge) = unkel(2,ie) - 2.0*(unkel(2,ie)*nx + &
                       unkel(3,ie)*ny)*nx        
         unkel(3,ge) = unkel(3,ie) - 2.0*(unkel(2,ie)*nx + &
                       unkel(3,ie)*ny)*ny        
         unkel(4,ge) = unkel(4,ie)

      endif

      if(case1 == 1) then           !Airfoil

         if(bface(3,ifac) == 4) then  !Far-field

            unkel(1,ge) = rhoin
            unkel(2,ge) = velin(1)
            unkel(3,ge) = velin(2)
            unkel(4,ge) = prein

         endif

     else if(case1 == 2) then     
        
        if(bface(3,ifac) == 4) then  !Inflow/Outflow

           c1 = sqrt(gam*unkel(4,ie)/unkel(1,ie))
           c0 = sqrt(gam*prein/rhoin)
           vn1 =  unkel(2,ie)*nx + unkel(3,ie)*ny
           vt1 = -unkel(2,ie)*ny + unkel(3,ie)*nx
           vn0 =  velin(1)*nx + velin(2)*ny
           vt0 = -velin(1)*ny + velin(2)*nx
           s0  = prein/rhoin**gam
           s1  = unkel(4,ie)/unkel(1,ie)**gam

           mn1 = vn1/c1

           if(mn1 < -1.0) then       !Supersonic inflow

              unkel(1,ge) = rhoin
              unkel(2,ge) = velin(1)
              unkel(3,ge) = velin(2)
              unkel(4,ge) = prein
           
           elseif (mn1 > 1.0) then   !Supersonic outflow    
            
              
              unkel(1,ge) = unkel(1,ie)
              unkel(2,ge) = unkel(2,ie)
              unkel(3,ge) = unkel(3,ie)
              unkel(4,ge) = unkel(4,ie)

           else if(mn1 >= -1.0 .and. mn1 < 0.0) then  !Subsonic inflow
   
              vng = 0.5*(vn0 + vn1) + (c1 - c0)/(gam - 1.0)
              vtg = vt0
              cg = 0.25*(gam - 1.0)*(vn1 -  vn0) + 0.5*(c1 + c0)

              unkel(1,ge) = rhoin*(cg/c0)**(2.0/(gam - 1.0))
              unkel(2,ge) = velin(1) - vn0*nx + vng*nx 
              unkel(3,ge) = velin(2) - vn0*ny + vng*ny 
              unkel(4,ge) = cg*cg*unkel(1,ge)/gam
        
           else if(mn1 > 0.0 .and. mn1 <= 1.0) then  !Subsonic outflow
   
              vng = 0.5*(vn0 + vn1) + (c1 - c0)/(gam - 1.0)
              vtg = vt1
              cg = 0.25*(gam - 1.0)*(vn1 -  vn0) + 0.5*(c0 + c1)

              unkel(1,ge) = unkel(1,ie)*(cg/c1)**(2.0/(gam - 1.0))
              unkel(2,ge) = unkel(2,ie) - vn1*nx + vng*nx 
              unkel(3,ge) = unkel(3,ie) - vn1*ny + vng*ny 
              unkel(4,ge) = cg*cg*unkel(1,ge)/gam
        
           endif

       endif

     endif

  enddo

end subroutine setbc
