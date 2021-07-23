!Module containing all the global variables
module comvar
   implicit none

   integer  ::  ndimn
   integer  ::  npoin, nelem, nafac, nbfac
   integer  ::  ndegr, neqns
   integer  ::  converge
   integer  ::  nnode, nconi
   integer  ::  nesup, mafac, mbfac
   integer  ::  nfael   !No.of faces per element
   integer  ::  case1

   real     ::  time
   real, allocatable, dimension(:,:)     ::  coord
   integer, allocatable, dimension(:,:)  ::  inpoel
   integer, allocatable, dimension(:,:)  ::  intfac
   integer, allocatable, dimension(:,:)  ::  bface
   integer, allocatable, dimension(:)    ::  esup1, esup2
   integer, allocatable, dimension(:,:)  ::  esuel
   
   real, allocatable, dimension(:,:)   ::  unkel !Primitive variables 
   real, allocatable, dimension(:,:,:)   ::  dunkel !Gradients of primitive variables
   real, allocatable, dimension(:,:)   ::  rhs 
   real, allocatable, dimension(:,:)   ::  cunkel !Conservative variables 
   real, allocatable, dimension(:,:)   ::  q0      
   real, allocatable, dimension(:,:)   ::  qn      
  
   real,allocatable, dimension(:,:)    ::  geoel
   real,allocatable, dimension(:,:)    ::  geofa
   
   real, allocatable, dimension(:,:)   ::  wmat
   real, allocatable, dimension(:) :: wgt, dt1

   real, parameter     ::  c_pi = 4.0*atan(1.0) 
   
   real  :: gam, gasR
   real  :: dt, cfl
   real  :: restol, ark(3)
   real   ::  rhoin, velin(2), prein
   integer :: itmax, itsave

   integer :: limtype, fluxtype
   real    :: lim_par ! Limiter value
   ! lim_par = 1  : Fromm's Method
   ! lim_par = -1 : Beam-Warming Method
   ! lim_par = 0  : Lax-Wendroff Method

   character(len = 100)   :: gfilename
   
end module comvar
