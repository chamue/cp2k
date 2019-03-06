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
                                               cp_fm_set_element, &
                                               cp_fm_to_fm, &
                                               cp_fm_create, &
                                               cp_fm_get_submatrix, &
                                               cp_fm_release
   USE cp_fm_struct,                     ONLY: cp_fm_struct_type, &
                                               cp_fm_struct_create, &
                                               cp_fm_struct_release
   USE lapack,                           ONLY: lapack_dgeqp3
   USE kinds,                            ONLY: dp
   USE cp_fm_basic_linalg,               ONLY: cp_fm_gemm                                  
   USE cp_para_types,                    ONLY: cp_para_env_type    
   USE message_passing,                  ONLY: mp_bcast   
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
  SUBROUTINE scdm_rs(vectors, para_env, sample)
   ! Input Variables
   TYPE(cp_fm_type), POINTER                     :: vectors
   TYPE(cp_para_env_type), POINTER               :: para_env
   logical, intent(in)                           :: sample
   !real(dp), intent(in), optional                :: eps

   ! Working Variables
#if defined(__SCALAPACK)
   type(cp_fm_type), pointer                     :: CTp, Qf, tmp
   INTEGER, DIMENSION(9)                         :: descc
   integer                                       :: lwork, ndim
   type(cp_fm_struct_type), pointer              :: cstruct
   real(kind=dp), dimension(:,:), pointer        :: a
   real(kind=dp), allocatable,  dimension(:,:)   :: matrix
   integer, allocatable, dimension(:)            :: ipiv
#endif
   integer                                       :: i, j, k, INFO, M, N, handle !, lld
   real(kind=dp)                                 :: norm
   real(kind=dp), allocatable, dimension(:,:)    :: Q, H, CTs, normmat !, v, Vtensor
   real(kind=dp), allocatable, dimension(:)      :: Tau, Work, roh, rohsum, U
   integer, allocatable, dimension(:)            :: jpvt 
   integer, allocatable, dimension(:,:)          :: C
   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN

   CALL timeset(routineN, handle)
   if (para_env%num_pe==1) then
        N = size(vectors%local_data,dim=2)
        M = size(vectors%local_data,dim=1)
        !print *, "N: ", N, ", M: ", M
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
        allocate(Q(N,N))
        allocate(H(N,N))
        allocate(jpvt(M))
  
        ! Pivoted QR Decomposition of Psicon
        jpvt = 0
        call lapack_dgeqp3(N,M,CTs,N,jpvt,Tau,Work,3*M+1,INFO)
        CPASSERT(INFO==0)

        ! Construction of Q from the lapack output
        do i = 1,N
          do j = 1,N
            do k=j,N
               if (j==k .and. j<i) then
                    H(j,k)=1.0_dp
               else if (j==k .and. j==i) then
                    H(j,k)=1.0_dp - tau(i)
               else if (j==k .and. j>i) then
                    H(j,k)=1.0_dp - tau(i)*CTs(j,i)*CTs(j,i)
               else if (j/=k .and. j==i .and. k>i) then
                    H(j,k)=-tau(i)*CTs(k,i)
                    H(k,j)=H(j,k)
               else if (j/=k .and. j>i .and. k>i) then
                    H(j,k)=-tau(i)*CTs(k,i)*CTs(j,i)
                    H(k,j)=H(j,k)
               else
                    H(k,j)=0.0_dp
                    H(j,k)=0.0_dp  
               end if 
            end do
          end do
          if (i == 1) then
            Q = H
          else
            Q = matmul(Q,H)
          end if
        end do
        deallocate(jpvt)
        ! Calculate localized orbitals
        vectors%local_data = matmul(real(vectors%local_data),real(Q))
        
        ! Clean up
        deallocate(CTs)
        deallocate(Tau)
        deallocate(Work)
        deallocate(Q)
        deallocate(H)

   else
   ! parallel version, no sampling 
#if defined(__SCALAPACK)
        if (sample) then 
           CPABORT("Sampling only defined for serial version")
        end if
        ! Create Transpose of C Matrix
        lwork=-1
        N = vectors%matrix_struct%ncol_global
        M = vectors%matrix_struct%nrow_global
        call cp_fm_struct_create(cstruct, vectors%matrix_struct%para_env, &
                vectors%matrix_struct%context, N, M, vectors%matrix_struct%ncol_block, &
                vectors%matrix_struct%nrow_block, first_p_pos=vectors%matrix_struct%first_p_pos)
        descc(:)=cstruct%descriptor(:)
        call cp_fm_create(CTp,cstruct)
        call pdtran(N,M,1.0_dp,vectors%local_data(1,1),1,1,vectors%matrix_struct%descriptor,&
                0.0_dp,CTp%local_data(1,1),1,1,descc)
        
        ! Allocating work variables
        allocate(Q(N,N))
        allocate(H(N,N))
        a => CTp%local_data
        ndim = size(a,2)
        allocate(tau(N))
        lwork=-1
        allocate(work(2*ndim))
        allocate(ipiv(M))

        ! Call SCALAPACK routine
        call pdgeqpf(n,m,a,1,1,descc,ipiv,tau,work,lwork,info)
        lwork=INT(work(1))
        deallocate(work)
        allocate(work(lwork))
        tau=0
        ipiv=0
        call pdgeqpf(n,m,a,1,1,descc,ipiv, tau,work,lwork,info)
        CPASSERT(INFO==0)
        ! Construction of Q from the scalapack output (should also be parallelized
        ! in the future)
        if (para_env%mepos==0 .and. all(tau==0)) then
          CPABORT("tau 0 for process 0")
        end if
        call mp_bcast(msg=tau,source=0,gid=para_env%group)
        allocate(matrix(N,N))
        call cp_fm_get_submatrix(CTp, matrix,1,1,N,N)
        do i = 1,N
          do j=1,N
            do k=j,N
               if (j==k .and. j<i) then
                    H(j,k)=1.0_dp
               else if (j==k .and. j==i) then
                    H(j,k)=1.0_dp - tau(i)
               else if (j==k .and. j>i) then
                    H(j,k)=1.0_dp - tau(i)*matrix(j,i)*matrix(j,i)
               else if (j/=k .and. j==i .and. k>i) then
                    H(j,k)=-tau(i)*matrix(k,i)
                    H(k,j)=H(j,k)
               else if (j/=k .and. j>i .and. k>i) then
                    H(j,k)=-tau(i)*matrix(k,i)*matrix(j,i)
                    H(k,j)=H(j,k)
               else
                    H(k,j)=0.0_dp
                    H(j,k)=0.0_dp
               end if
            end do
          end do
          if (i == 1) then
            Q = H
          else
            Q=matmul(Q,H)
          end if   
        end do
        
        ! Construct fm Q
        call cp_fm_struct_release(cstruct)
        call cp_fm_struct_create(cstruct,CTp%matrix_struct%para_env, CTp%matrix_struct%context,&
               CTp%matrix_struct%nrow_global, CTp%matrix_struct%nrow_global, &
               CTp%matrix_struct%nrow_block, CTp%matrix_struct%nrow_block, &
               local_leading_dimension=CTp%matrix_struct%local_leading_dimension)
        call cp_fm_create(Qf,cstruct)
        do i = 1,N
          do j = 1,N
             call cp_fm_set_element(Qf,i,j,Q(i,j))
          end do
        end do
        
        ! Transform original coefficient matrix vectors
        call cp_fm_create(tmp,vectors%matrix_struct)
        call cp_fm_to_fm(vectors, tmp)
        call cp_fm_gemm('N','N',vectors%matrix_struct%nrow_global,vectors%matrix_struct%ncol_global,&
                vectors%matrix_struct%ncol_global,1.0_dp,tmp,Qf,0.0_dp,vectors)
        
        ! Clean up
        call cp_fm_struct_release(cstruct)
        call cp_fm_release(CTp)
        call cp_fm_release(Qf)
        deallocate(tau)
        deallocate(work)
        deallocate(Q)
        deallocate(H)
        deallocate(ipiv)
        deallocate(matrix)

#else
        CPABORT("pdgeqrf() called without SCALAPACK compiled in")
#endif
   end if 
   CALL timestop(handle)
   
  END SUBROUTINE scdm_rs
end module scdm
