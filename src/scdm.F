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
  SUBROUTINE scdm_rs(vectors)
   ! Input Variables
   TYPE(cp_fm_type), POINTER                     :: vectors

   ! Working Variables
   integer                                       :: i, j, k, INFO, M, N, handle
   real(kind=dp), allocatable, dimension(:,:)    :: Q, H, Vtensor,CT
   real(kind=dp), allocatable, dimension(:)      :: Tau, v, Work
   integer, allocatable, dimension(:)            :: jpvt 
   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN

   CALL timeset(routineN, handle)
   
   ! Allocating work variables
   N = size(vectors%local_data,dim=2)
   M = size(vectors%local_data,dim=1)
   allocate(Tau(N))
   allocate(v(N))
   allocate(Work(3*M+1))
   allocate(Q(N,N))
   allocate(H(N,N))
   allocate(jpvt(M))
   allocate(Vtensor(N,N))
   allocate(CT(N,M))
   ! Pivoted QR Decomposition of Psicon
   jpvt = 0
   CT = transpose(vectors%local_data)
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

