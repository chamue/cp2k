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
   USE cp_fm_types,                      ONLY: cp_fm_type,&
                                               cp_fm_to_fm, &
                                               cp_fm_create, &
                                               cp_fm_release, &
                                               cp_fm_to_fm_submat
   USE cp_fm_struct,                     ONLY: cp_fm_struct_type, &
                                               cp_fm_struct_create, &
                                               cp_fm_struct_release
   USE lapack,                           ONLY: lapack_dgeqp3
   USE kinds,                            ONLY: dp
   USE cp_fm_basic_linalg,               ONLY: cp_fm_gemm                                  
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
  SUBROUTINE scdm_rs(vectors, sample)
   ! Input Variables
   TYPE(cp_fm_type), POINTER                     :: vectors
   logical, intent(in)                           :: sample

   ! Working Variables
#if defined(__SCALAPACK)
   integer                                       :: INFO, M, N, handle
   real(kind=dp), allocatable, dimension(:)      :: Tau, Work
   type(cp_fm_type), pointer                     :: CTp, tmp, Qf
   INTEGER, DIMENSION(9)                         :: descc
   integer                                       :: lwork
   type(cp_fm_struct_type), pointer              :: cstruct
   real(kind=dp), dimension(:,:), pointer        :: a
   integer, allocatable, dimension(:)            :: ipiv
#else
   integer                                       :: i, j, lwork, INFO, M, N, handle
   real(kind=dp)                                 :: norm
   real(kind=dp), allocatable, dimension(:,:)    :: CTs, normmat 
   real(kind=dp), allocatable, dimension(:)      :: Tau, Work, roh, rohsum, U
   integer, allocatable, dimension(:)            :: jpvt 
   integer, allocatable, dimension(:,:)          :: C
#endif
   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN
   CALL timeset(routineN, handle)
#if defined(__SCALAPACK)
        if (sample) then
           CPABORT("Sampling only defined for serial version")
        end if
        ! Create Transpose of Coefficient Matrix vectors
        N = vectors%matrix_struct%ncol_global
        M = vectors%matrix_struct%nrow_global
        call cp_fm_struct_create(cstruct, vectors%matrix_struct%para_env, &
                vectors%matrix_struct%context, N, M, vectors%matrix_struct%ncol_block, &
                vectors%matrix_struct%nrow_block, first_p_pos=vectors%matrix_struct%first_p_pos)
        descc(:)=cstruct%descriptor(:)
        call cp_fm_create(CTp,cstruct)
        call pdtran(N,M,1.0_dp,vectors%local_data(1,1),1,1,vectors%matrix_struct%descriptor,&
                0.0_dp,CTp%local_data(1,1),1,1,descc)

        ! Prepare work variables
        a => CTp%local_data
        allocate(tau(N))
        lwork=-1
        allocate(work(2*N))
        allocate(ipiv(M))
        
        ! Call SCALAPACK routine to get optimal work dimension
        call pdgeqpf(n,m,a,1,1,descc,ipiv,tau,work,lwork,info)
        lwork=INT(work(1))
        deallocate(work)
        allocate(work(lwork))
        tau=0.0_dp
        ipiv=0
        ! Call SCALAPACK routine to get QR decomposition of CTs
        call pdgeqpf(n,m,a,1,1,descc,ipiv, tau,work,lwork,info)
        CPASSERT(INFO==0)

        ! Construction of Q from the scalapack output
        call cp_fm_struct_release(cstruct)
        call cp_fm_struct_create(cstruct,para_env=CTp%matrix_struct%para_env, &
                context=CTp%matrix_struct%context, nrow_global=CTp%matrix_struct%nrow_global, &
                ncol_global=CTp%matrix_struct%nrow_global)
        call cp_fm_create(Qf,cstruct)
        call cp_fm_to_fm_submat(CTp,Qf,N,N,1,1,1,1)
        a => Qf%local_data
        lwork=-1
        ! Call SCALAPACK routine to get optimal work dimension
        call pdorgqr(n,n,n,a,1,1,descc,tau,work,lwork,info)
        CPASSERT(info==0)
        lwork=int(work(1))
        deallocate(work)
        allocate(work(lwork))
        ! Call SCALAPACK routine to get Q 
        call pdorgqr(n,n,n,a,1,1,descc,tau,work,lwork,info)
        CPASSERT(info==0)

        ! Transform original coefficient matrix vectors
        call cp_fm_create(tmp,vectors%matrix_struct)
        call cp_fm_to_fm(vectors, tmp)
        call cp_fm_gemm('N','N',M,N,N,1.0_dp,tmp,Qf,0.0_dp,vectors)
        
        ! Clean up
        call cp_fm_struct_release(cstruct)
        call cp_fm_release(CTp)
        call cp_fm_release(Qf)
        deallocate(tau)
        deallocate(work)
        deallocate(ipiv)
#else
        N = size(vectors%local_data,dim=2)
        M = size(vectors%local_data,dim=1)
        ! This sampling is unnecessary when working with the coefficient matrix
        ! It was created for the case when working in realspace
        ! Maybe use the sample keyword in the future to see if BS should be 
        ! projected on a minimal basis when not well behaved
        if (sample) then
           allocate(roh(M))
           allocate(rohsum(M))
           norm = 0
           allocate(normmat(N,N))
           normmat=matmul(transpose(vectors%local_data),vectors%local_data)
           do i = 1,N
             norm=norm+normmat(i,i)
           end do
           ! Calculate probability distribution
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
           rohsum = rohsum/rohsum(M)
           ! Random Sample
           M=ceiling(3*N*log(real(N)))
           allocate(C(M,1))
           allocate(U(M))
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
           allocate(CTs(N,M))
           CTs = transpose(vectors%local_data(C(:,1),:))
           deallocate(U)
           deallocate(C)
           deallocate(normmat)
           deallocate(rohsum)
           deallocate(roh)
        else
           allocate(CTs(N,M))
           CTs = transpose(vectors%local_data)
        end if

        ! Allocating work variables
        allocate(Tau(N))
        allocate(Work(3*M+1))
        allocate(jpvt(M))

        ! Pivoted QR Decomposition of Psicon
        jpvt = 0
        call lapack_dgeqp3(N,M,CTs,N,jpvt,Tau,Work,3*M+1,INFO)
        CPASSERT(INFO==0)
        ! Get Q from QR output
        call dorgqr(N,N,N,CTs(1:N,1:N),N,Tau,work,-1,info)
        lwork=int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dorgqr(N,N,N,CTs(1:N,1:N),N,Tau,work,lwork,info)
        CPASSERT(info==0)
        
        ! Calculate localized orbitals
        vectors%local_data = matmul(real(vectors%local_data),real(CTs(1:N,1:N)))

        ! Clean up
        deallocate(CTs)
        deallocate(Tau)
        deallocate(Work)
        deallocate(jpvt)

#endif
   CALL timestop(handle)
  END SUBROUTINE scdm_rs
end module scdm
