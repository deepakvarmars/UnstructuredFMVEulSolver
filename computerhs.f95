subroutine computerhs()
   use comvar
   implicit none
   integer :: j, ifac, ielem, jelem
   real    :: lng, nx, ny
   real    :: ql(4), qr(4), flux(4)
   real    :: qi(4), qj(4), dqi(2,4), dqj(2,4), dr(2)
   
   rhs = 0.0d0
   dunkel = 0.0d0
   dqi = 0.0d0
   dqj = 0.0d0
   flux = 0.0d0

   do ifac = 1, nafac

      lng = geofa(3,ifac)             !length of the face/edge
      nx  = geofa(1,ifac)/lng         !Unit normal - x direction
      ny  = geofa(2,ifac)/lng         !Unit normal - y direction
  
      ielem = intfac(1,ifac)
      jelem = intfac(2,ifac)

      if(limtype > 1) then

         call gradcompute()

         if(ifac <= nbfac ) then  !boundary face

            dr(1) = geofa(1,ifac)
            dr(2) = geofa(2,ifac)

         else
 
            dr(1) = geoel(1,jelem) - geoel(1,ielem)
            dr(2) = geoel(2,jelem) - geoel(2,ielem)
         endif

      endif
      qi(:)    = unkel(:,ielem)
      qj(:)    = unkel(:,jelem)
      dqi(:,:) = dunkel(:,:,ielem)
      dqj(:,:) = dunkel(:,:,jelem)

      call reconstruct(qi, qj, dqi, dqj, dr, ql, qr)

      ql(4) = max(0.0000001d0, ql(4))
      ql(1) = max(0.0000001d0, ql(1))
      qr(4) = max(0.0000001d0, qr(4))
      qr(1) = max(0.0000001d0, qr(1))

      call numflux(nx, ny, lng, ql, qr, flux)

      do j=1,4
         if(isnan(flux(j))) then
            write(*,*) "flux NaN"
            write(*,*) "For face no ", ifac
            write(*,*) "bface type", bface(3,ifac)
            write(*,*) ielem, jelem
            call exit(0)
         endif
      enddo

      rhs(:,ielem) = rhs(:,ielem) - flux(:)
      rhs(:,jelem) = rhs(:,jelem) + flux(:)

   enddo

   return

end subroutine computerhs

subroutine reconstruct(qi,qj,dqi,dqj,dr,ql,qr)
    use comvar
    implicit none
    real, intent(in)  :: qi(4), qj(4), dr(2)
    real, intent(in)  :: dqi(2, 4), dqj(2, 4)
    real, intent(out) :: ql(4), qr(4)
    !Local variables
    real  :: eps = 1.0d-8
    real  :: dq(4), dtmpl(4), dtmpr(4)
    real  :: phil(4), phir(4)
    
    if(limtype == 1) then    !First order
       ql(:) = qi(:)
       qr(:) = qj(:)
    else if(limtype == 2 ) then     !Second Order MUSCL with van Albada
       
       dq(:) = qj(:) - qi(:)
       
       dtmpl(:) = 2.0d0*(dqi(1,:)*dr(1) + dqi(2,:)*dr(2)) - dq(:)
       dtmpr(:) = 2.0d0*(dqj(1,:)*dr(1) + dqj(2,:)*dr(2)) - dq(:)

       !van Albada limiter
       phil(:) = max(0.0d0, 2.0*dtmpl(:)*dq(:)/&
                 (dtmpl(:)**2 + dq(:)**2 + eps))
       phir(:) = max(0.0d0, 2.0*dtmpr(:)*dq(:)/&
                 (dtmpr(:)**2 + dq(:)**2 + eps))
       
       ql(:) = qi(:) + 0.25d0*phil(:)*((1.0d0 - lim_par*phil(:))*dtmpl(:) + &
               (1.0d0 + lim_par*phil(:))*dq(:))
       qr(:) = qj(:) - 0.25d0*phir(:)*((1.0d0 - lim_par*phir(:))*dtmpr(:) + &
               (1.0d0 + lim_par*phir(:))*dq(:))

    endif

    return
 
end subroutine reconstruct

subroutine numflux(nx, ny, lng, ql, qr, flux)         !l = left = i; r = right = j
    
   use comvar
   implicit none
   real, intent(in)  :: nx, ny, lng, ql(4), qr(4)
   real, intent(out) :: flux(4)
   !local variables 
   !integer :: i, j, k
   !Variables for van Leer flux
   real  :: vni, vnj
   real  :: flxl(4), flxr(4)
   real  :: mni, mnj, ci, cj 
   real  :: ri, ui, vi, pi, ei
   real  :: rj, uj, vj, pj, ej

!!!!INSERT FLUX SCHEME HERE
   
   vni = ql(2)*nx + ql(3)*ny
   vnj = qr(2)*nx + qr(3)*ny
   ci = sqrt(gam*ql(4)/ql(1))
   cj = sqrt(gam*qr(4)/qr(1))
   mni = vni/ci
   mnj = vnj/cj
   ri = ql(1) ; rj = qr(1);
   ui = ql(2) ; uj = qr(2);
   vi = ql(3) ; vj = qr(3);
   pi = ql(4) ; pj = qr(4);

   ei = pi/(ri*(gam-1.0)) + 0.5*(ui**2 + vi**2)
   ej = pj/(rj*(gam-1.0)) + 0.5*(uj**2 + vj**2)

   if (fluxtype == 2) then   !van Leer flux-splitting scheme

      ! calculation of left flux
      if ( mni > 1.0d0 ) then

        flxl(1) = ri*vni
        flxl(2) = ri*vni*ui + pi*nx
        flxl(3) = ri*vni*vi + pi*ny
        flxl(4) = vni*(ri*ei + pi)

      elseif ( mni < -1.0d0 ) then

        flxl = 0.0d0

      else

        flxl(1) =  0.25d0*ri*ci*(mni+1.0d0)**2
        flxl(2) = flxl(1)*( ui + nx*(-vni + 2.0d0*ci)/gam )
        flxl(3) = flxl(1)*( vi + ny*(-vni + 2.0d0*ci)/gam )
        flxl(4) = flxl(1)*( 0.5d0*(ui**2 + vi**2 - vni**2) &
          + 0.5d0*( ((gam-1.0)*vni + 2.0d0*ci)**2 )/(gam**2 - 1.0d0) )

      endif

      ! calculation of right flux
      if ( mnj > 1.0d0 ) then

        flxr = 0.0d0

      elseif ( mnj < -1.0d0 ) then

        flxr(1) = rj*vnj
        flxr(2) = rj*vnj*uj + pj*nx
        flxr(3) = rj*vnj*vj + pj*ny
        flxr(4) = vnj*(rj*ej + pj)

      else

        flxr(1) = -0.25d0*rj*cj*(mnj-1.0d0)**2
        flxr(2) = flxr(1)*( uj + nx*(-vnj - 2.0d0*cj)/gam )
        flxr(3) = flxr(1)*( vj + ny*(-vnj - 2.0d0*cj)/gam )
        flxr(4) = flxr(1)*( 0.5d0*(uj**2 + vj**2 - vnj**2) &
          + 0.5d0*( ((gam-1.0)*vnj - 2.0d0*cj)**2 )/(gam**2 - 1.0d0) )

      endif
      
      flux(:) = (flxl(:) + flxr(:))*lng
   
     endif

   return

end subroutine numflux

subroutine gradcompute()
   use comvar
   implicit none
   integer :: ielem, ifael, jelem, ifac
   real    :: coeff, wj, dxj, dyj
   real    :: matr(ndimn,4)

   do ielem = 1, nelem
  
      matr = 0.0
     
      do ifael = 1,nfael

         jelem = esuel(ifael,ielem)

         if(jelem > nelem) then

            ifac = jelem - nelem
            dxj = geofa(1,ifac)
            dyj = geofa(2,ifac)

            wj = 1.0d0/geofa(3,ifac)**2

         else

            dxj = geoel(1,jelem) - geoel(1,ielem)
            dyj = geoel(2,jelem) - geoel(2,ielem)

            wj = 1.0d0/(dxj**2 + dyj**2)

         endif

         matr(1,:) = matr(1,:) + wj*dxj*(unkel(:,jelem) - unkel(:,ielem))
         matr(2,:) = matr(2,:) + wj*dyj*(unkel(:,jelem) - unkel(:,ielem))
     
      enddo

      !Obtaining the grads
      coeff = (wmat(1,ielem)*wmat(3,ielem) - wmat(2,ielem))**2    !Determinant of LHS matrix
 
       dunkel(1,:,ielem) = (matr(1,:)*wmat(3,ielem) - matr(2,:)*wmat(2,ielem))/coeff
       dunkel(2,:,ielem) = (matr(1,:)*wmat(2,ielem) - matr(2,:)*wmat(1,ielem))/coeff

   enddo

   !For BCs, ghost cell extrapolation is required

    do ifac = 1, nbfac

       ielem = intfac(1,ifac)
       jelem = intfac(2,ifac)


       dunkel(:,:,jelem) = dunkel(:,:,ielem)

    enddo 
   
    return 
 
end subroutine gradcompute
