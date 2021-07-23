subroutine writesol(fnum)
   use comvar
   implicit none
   integer, intent(in)  :: fnum
   integer  ::  ipoin, ielem, icoun, i
   real     :: wgt1, m, velmag, ss1
   real     :: s_in, sn, err1, gridval
   character(len=50)  :: filenm

   write(*,*) "Postprocessing and writing solution file"
  
   qn = 0.0

   do ipoin = 1, npoin

      wgt1 = 0.0

      do icoun = esup2(ipoin) + 1, esup2(ipoin+1)

        ielem = esup1(icoun)

        !Area weighted sum of the elements surrounding the point
        qn(:,ipoin) = qn(:,ipoin) + unkel(:,ielem)*geoel(3,ielem)
        
        wgt1 = wgt1 + geoel(3,ielem)   !Total area of surrounding elements

      enddo

      qn(:,ipoin) = qn(:,ipoin)/wgt1

   enddo

   filenm = "Solution"
   write(unit=filenm,fmt='(A8,I0.5)') trim(filenm), fnum
   filenm = trim(filenm)//'.plt'

   open(unit = 50, file = trim(filenm))
   write(50,*) 'Variables = "X","Y","Rho","U","V","Pre","M"'
   write(50,*) "Zone F = fepoint, N =", npoin,", E = ", nelem,", ET = Triangle"

   do ipoin = 1, npoin
      ss1 = sqrt(gam*qn(4,ipoin)/qn(1,ipoin))
      velmag = sqrt(qn(2,ipoin)**2 + qn(3,ipoin)**2)
      m = velmag/ss1
      write(50,'(7e30.12)') (coord(i,ipoin),i=1,2), qn(1,ipoin), qn(2,ipoin)&
                            , qn(3,ipoin), qn(4,ipoin), m
   enddo

   do ielem=1, nelem
      write(50,*) (inpoel(i,ielem), i=1,3)
   enddo

   err1 = 0.0d0
   s_in = prein/rhoin**gam

   do ipoin=1,npoin

      sn = qn(4,ipoin)/qn(1,ipoin)**gam

      err1 = err1 + ((sn - s_in)/s_in)**2

   enddo

   err1 = sqrt(err1/npoin)

   gridval = 0.0d0

   do ielem = 1, nelem
      gridval = gridval + geoel(3,ielem)
   enddo

   gridval = sqrt(gridval/nelem)

   write(*,*) "Error, Grid_Value"
   write(*,*) err1, gridval
   close(50)


end subroutine writesol
