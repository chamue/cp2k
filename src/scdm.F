!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Localization method scdm using pivoted QR decomposition
!> \par History Implementation scdm during master thesis
!>              06.12.18: Serial Implementation Done, Commented out parts are 
!>                        mostly for parallel uses in the future
!>
!> \author Charlotte Mueller (12.2018) 
! *****************************************************************************
MODULE scdm
   USE cp_fm_types,                      ONLY: cp_fm_type
   USE lapack,                           ONLY: lapack_dgeqp3
   USE kinds,                            ONLY: dp

#include "./base/base_uses.f90"
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: scdm_rs

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scdm'

CONTAINS
! *****************************************************************************
!> \brief Computes localized orbitals by pivoted QR decomposition. Might be 
!>        optimized by random sampling.
!> \param vectors: containing non-localized mo-coefficients
!> \param weights: root of volume of single voxel (in realspace)
!>
! *****************************************************************************
  SUBROUTINE scdm_rs(vectors, sample, eps)
   ! Input Variables
   TYPE(cp_fm_type), POINTER                     :: vectors
   logical, intent(in)                           :: sample
   real(dp), intent(in), optional                :: eps

   ! Working Variables
   integer                                       :: i, j, k, INFO, M, N, handle
   real(kind=dp), allocatable, dimension(:,:)    :: Q, H, Vtensor,CT
   real(kind=dp), allocatable, dimension(:)      :: Tau, v, Work, roh, rohsum
   integer, allocatable, dimension(:)            :: jpvt 
   integer, dimension(100,1)                     :: C
   real, dimension(100)                          :: U
   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN

   CALL timeset(routineN, handle)
   N = size(vectors%local_data,dim=2)
   M = size(vectors%local_data,dim=1)
   allocate(roh(M))
   allocate(rohsum(M))
   if (sample) then
     ! Calculate probability distribution
     print *, "Here"
     do i = 1,M
       roh(i) = 0
       do j = 1,N
         roh(i) = roh(i)+vectors%local_data(i,j)*vectors%local_data(i,j)
       end do
       if (i == 1) then
         rohsum(i) = roh(i)
       else 
         rohsum(i) = rohsum(i-1) + roh(i)
       end if 
     end do
     print *, rohsum(M)
     rohsum = rohsum/rohsum(M)
     ! Random Sample
     call random_number(U)
     C = 0
     do i = 1,size(U)
       do j = 1, size(roh)
         if (rohsum(j-1) <= U(i) .AND. U(i) <= rohsum(j) .AND. C(i,1) == 0) then
           C(i,1) = j
           exit
         end if
       end do
     end do
     M = 100
     allocate(CT(N,M))
     CT = transpose(vectors%local_data(C(:,1),:))
   else 
     allocate(CT(N,M))
     CT = transpose(vectors%local_data)  
   end if 
   ! Allocating work variables
   allocate(Tau(N))
   allocate(v(N))
   allocate(Work(3*M+1))
   allocate(Q(N,N))
   allocate(H(N,N))
   allocate(jpvt(M))
   allocate(Vtensor(N,N))
  
   ! Pivoted QR Decomposition of Psicon
   jpvt = 0
   call lapack_dgeqp3(N,M,CT,N,jpvt,Tau,Work,3*M+1,INFO)
   CPASSERT(INFO==0)

   ! Construction of Q from the lapack output
   do i = 1,N
     H = 0
     v(1:i-1)=0
     v(i)=1
     v(i+1:N)=CT(i+1:N,i)
     do j = 1,N
       H(j,j)=1
       do k = 1,N
         Vtensor(j,k) = v(j)*v(k)
       end do
     end do
     H = H - Tau(i)*Vtensor
     if (i == 1) then
       Q = H
     else
       !print *, "Q*H: Q:", shape(Q), "H:", shape(H)
       Q = matmul(Q,H)
     end if
   end do
   
   ! Calculate localized orbitals
   !print *, "local_data*Q: local_data:", shape(real(vectors%local_data)), &
   !        "Q:", shape(Q)
   vectors%local_data = matmul(real(vectors%local_data),real(Q))
   
   CALL timestop(handle)
   
  END SUBROUTINE scdm_rs

end module scdm

