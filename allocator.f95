subroutine allocator
   use comvar
   implicit none

   mafac = 3*nelem
   nfael = 3
   neqns = 4

   allocate(coord(ndimn, npoin))
   allocate(inpoel(nnode, nelem))
   allocate(bface(4, nbfac))

   allocate(esup2(npoin+1))
   allocate(esuel(nnode, nelem))

   allocate(unkel(4, nelem+nbfac))    !Primitive variables
   allocate(cunkel(4, nelem+nbfac))   !Conservative  variables
   allocate(dunkel(2,4, nelem+nbfac))   !Gradients of primitive variables
   allocate(rhs(4, nelem+nbfac))
   allocate(q0(4, nelem+nbfac))    
   allocate(qn(4, npoin))    

   allocate(wmat(3,nelem))
   
   allocate(wgt(nelem), dt1(nelem))

end subroutine allocator
