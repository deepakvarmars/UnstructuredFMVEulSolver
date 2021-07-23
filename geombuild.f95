subroutine geombuild
  use comvar
  implicit none

  integer   :: ielem,ifac
  integer   :: ip1,ip2,ip3
  real(8)   :: x1,x2,x3,y1,y2,y3,a1,b1,a2,b2
  real(8)   :: nx,ny,lng
  real(8), parameter :: obt = 1.0d0/3.0d0

  allocate(geoel(3,nelem))
  allocate(geofa(3,nafac))

  do ielem = 1,nelem

    ip1 = inpoel(1,ielem)
    ip2 = inpoel(2,ielem)
    ip3 = inpoel(3,ielem)

    x1 = coord(1,ip1); y1 = coord(2,ip1);
    x2 = coord(1,ip2); y2 = coord(2,ip2);
    x3 = coord(1,ip3); y3 = coord(2,ip3);

    ! Compute  centroid
    geoel(1,ielem) = obt*(x1 + x2 + x3)
    geoel(2,ielem) = obt*(y1 + y2 + y3)

    a1 =   y2 - y3
    b1 = -(x2 - x3)

    a2 =   y3 - y1
    b2 = -(x3 - x1)

    geoel(3,ielem) = 0.5*( a1*b2 - a2*b1 ) ! volume of the element

  enddo

  do ifac = 1,nafac

      ip1 = intfac(3,ifac)
      ip2 = intfac(4,ifac)

      x1 = coord(1,ip1); y1 = coord(2,ip1);
      x2 = coord(1,ip2); y2 = coord(2,ip2);

      nx = y2 - y1  ! normals
      ny = x1 - x2

      lng = sqrt(nx**2 + ny**2)  ! length of the face

      ! geometry of the faces
      geofa(1,ifac) = nx
      geofa(2,ifac) = ny
      geofa(3,ifac) = lng

  enddo

end subroutine geombuild
