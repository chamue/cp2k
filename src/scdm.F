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
   USE atomic_kind_types,                ONLY: get_atomic_kind
   USE qs_environment_types,             ONLY: get_qs_env, &
                                               qs_environment_type
   USE atomic_kind_types,                ONLY: atomic_kind_type
   USE qs_kind_types,                    ONLY: qs_kind_type
   USE cell_types,                       ONLY: cell_type
   USE cp_control_types,                 ONLY: dft_control_type
   USE qs_mo_types,                      ONLY: mo_set_p_type
   USE particle_types,                   ONLY: particle_type
   USE particle_list_types,              ONLY: particle_list_type
   USE qs_subsys_types,                  ONLY: qs_subsys_get,&
                                               qs_subsys_type
   USE pw_env_types,                     ONLY: pw_env_get,&
                                               pw_env_type    
   USE pw_pool_types,                    ONLY: pw_pool_type, &
                                               pw_pool_p_type, &
                                               pw_pool_create_pw, &
                                               pw_pool_give_back_pw
   USE pw_types,                         ONLY: COMPLEXDATA1D,&
                                               REALDATA3D,&
                                               REALSPACE,&
                                               RECIPROCALSPACE,&
                                               pw_p_type, &
                                               pw_type
   USE qs_mo_types,                      ONLY: get_mo_set
   USE cp_fm_types,                      ONLY: cp_fm_type, &
                                               cp_fm_create, &
                                               cp_fm_set_all, &
                                               cp_fm_release, &
                                               cp_fm_p_type
   USE qs_collocate_density,             ONLY: calculate_wavefunction
   !USE message_passing,                  ONLY: mp_maxloc,&
   !                                            mp_recv,&
   !                                            mp_send,&
   !                                            mp_sum, &
   !                                            mp_sync
   USE lapack,                          ONLY: lapack_dgeqp3
   USE kinds,                           ONLY: dp
   USE qs_localization_methods,         ONLY: rotate_orbitals
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
  
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
  SUBROUTINE scdm_rs(qs_env, zij, vectors, ispin, stride)
   ! Input Variables
   TYPE(qs_environment_type), POINTER            :: qs_env
   TYPE(cp_fm_p_type), INTENT(INOUT)             :: zij(:,:)
   TYPE(cp_fm_type), POINTER                     :: vectors
   INTEGER                                       :: ispin
   INTEGER, DIMENSION(:), OPTIONAL, POINTER      :: stride

   ! Working Variables
   integer                                       :: k, i, j, INFO, reason, M, N, &
                                                    ivector, ne, np
   real(kind=dp), allocatable, dimension(:,:)    ::  vectors_rs_2d ! M, N
   real(kind=dp), allocatable, dimension(:,:,:,:) :: vectors_rs_4d
   double precision, allocatable, dimension(:)   :: Tau, v ! N
   double precision, allocatable, dimension(:)   :: Work ! 3M+1
   double precision, allocatable, dimension(:,:) :: Q, H, Vtensor ! N, N
   double precision                              :: omega
   integer, allocatable, dimension(:)            :: jpvt 
   INTEGER :: homo, my_stride
   TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                    :: atomic_kind_set !!! get_qs_env
   TYPE(cp_fm_type), POINTER                     :: mo_coeff        !!! get_mo_set
   TYPE(dft_control_type), POINTER               :: dft_control     !!! get_qs_env
   TYPE(mo_set_p_type), DIMENSION(:), &
       POINTER                                   :: mos             !!! get_qs_env
   TYPE(particle_list_type), POINTER             :: particles       !!! particle_list_create
   TYPE(particle_type), DIMENSION(:), &
      POINTER                                    :: particle_set    !!! get_qs_env
   TYPE(pw_env_type), POINTER                    :: pw_env          !!! get_qs_env
   TYPE(pw_p_type)                               :: wf_g, wf_r      
   TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                    :: qs_kind_set     !!! get_qs_env
   TYPE(particle_list_type), POINTER             :: my_particles    !!! my_particles => particles
   TYPE(pw_pool_type), POINTER                   :: auxbas_pw_pool  !!! pw_env_get
   TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                    :: pw_pools        !!! pw_env_get
   TYPE(cell_type), POINTER                      :: cell            !!! get_qs_env
   INTEGER :: checksum, dest, gid, I1, I2, I3, N2, N3, ind, ip, &
     my_rank, num_pe, rank(2), source, stat, tag
   LOGICAL                                       :: failure
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)      :: buf
   type(cp_fm_type), pointer                     :: rmat           
   TYPE(qs_subsys_type), POINTER                 :: subsys
   integer, dimension(3)                         :: L, U, Lmax, Umax

   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN

   ! Get working variables
   NULLIFY(dft_control, pw_env, auxbas_pw_pool, pw_pools, mos, mo_coeff, &
           subsys, particles, my_particles)
   CALL get_qs_env(qs_env=qs_env,&
                         atomic_kind_set=atomic_kind_set,&
                         qs_kind_set=qs_kind_set,&
                         cell=cell,&
                         dft_control=dft_control,&
                         mos=mos, &
                         subsys=subsys, &
                         particle_set=particle_set, &
                         pw_env=pw_env)
   CALL qs_subsys_get(subsys, particles=particles)
   CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                          pw_pools=pw_pools)
   CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE)
   ! When no wf_g, error 'Missing actual argument' in calculate_wavefunction
   CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE)
   CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
   my_particles=>particles
   ne = my_particles%n_els
   failure=.FALSE.
   my_stride = 2
   
   ! for now, scdm_rs is called without stride input variable
   ! This code part is for future use of stride input
   if (present(stride)) then 
      CPASSERT(SIZE(stride)==1.OR.SIZE(stride)==3)
      IF (SIZE(stride)==3) THEN
         CPASSERT(stride(1)==stride(2).AND.stride(1)==stride(3))
      END IF  
      my_stride = stride(1)
   end if
   CPASSERT(my_stride>0)
   ! Looping over occupied orbitals and inserting real-space values in matrix columns
   DO ivector=1,homo
         CALL calculate_wavefunction(mo_vectors=mo_coeff,ivector=ivector,rho=wf_r,&
                 rho_gspace=wf_g,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
              cell=cell,dft_control=dft_control,particle_set=particle_set,&
              pw_env=pw_env)
          omega = wf_r%pw%pw_grid%dh(1,1)*REAL(my_stride,dp)
          omega = omega*omega*omega
          omega = sqrt(omega)
          do i = 1,3
            L(i) = wf_r%pw%pw_grid%bounds(1,i)
            U(i) = wf_r%pw%pw_grid%bounds(2,i)
          end do
          IF (ivector==1) then
            !allocate(vectors_rs_4d(L(3):U(3),L(2):U(2),L(1):U(1),homo)) 
            !ALLOCATE(vectors_rs_2d(homo,(U(1)-L(1)+1)*(U(2)-L(2)+1)*(U(3)-L(3)+1)))
            allocate(vectors_rs_4d((U(1)-L(1))/2+1,(U(2)-L(2))/2+1,(U(3)-L(3))/2+1,homo))
            ALLOCATE(vectors_rs_2d(homo,((U(3)-L(3))/2+1)*((U(2)-L(2))/2+1)*((U(1)-L(1))/2+1)))
          END IF
           do i=L(3),U(3),2
             do j=L(2),U(2),2
               do k=L(1),U(1),2
                 ! 
                 vectors_rs_4d((k-L(1))/2+1,(j-L(2))/2+1,(i-L(3))/2+1,ivector)=wf_r%pw%cr3d(k,j,i)
               end do
             end do
           end do
           do k=1,size(vectors_rs_4d,dim=3)
             do j = 1,size(vectors_rs_4d,dim=2)
               do i = 1,size(vectors_rs_4d, dim=4)
                  if (vectors_rs_4d(i,j,k,ivector)==0) print *, i,j,k,ivector
               end do
             end do
           end do
           !vectors_rs_4d(ivector,:,:,:)=wf_r%pw%cr3d(:,:,:)
   END DO

   ! Normalize vectors_rs
   vectors_rs_4d = vectors_rs_4d
   vectors_rs_2d=transpose(reshape(vectors_rs_4d, (/ size(vectors_rs_4d)/homo,homo /)))
   print *, "vectors_rs_2d: ", shape(vectors_rs_2d)
   print *, vectors_rs_2d(1,:)
   print *, "start scdm calculation"
   print *, "vectors_rs_4d: ", shape(vectors_rs_4d)
   ! Allocating work variables
   M = size(vectors_rs_2d,dim=2)
   N = size(vectors_rs_2d,dim=1)
   allocate(Tau(N))
   allocate(v(N))
   allocate(Work(3*M+1))
   allocate(Q(N,N))
   allocate(H(N,N))
   ! Pivoted QR Decomposition of Psicon
   call lapack_dgeqp3(N,M,vectors_rs_2d,N,jpvt,Tau,Work,3*M+1,INFO)
   CPASSERT(INFO==0)
   do i = 1,N
      H = 0
      v(1:i-1)=0
      v(i)=1
      v(i+1:N)=vectors_rs_2d(i+1:N,i)
      do j = 1,N
         H(j,j)=1
      end do
      do k = 1,N
        do j = 1,N
          Vtensor(k,j) = v(k)*v(j)
        end do
      end do
      H = H - Tau(i)*Vtensor
      if (i == 1) then
         Q = H
      else
         Q = matmul(Q,H)
      end if
   end do
!   Q =reshape((/-0.31355579783702092,&     
!   3.3047567236581368E-008,&
! -0.21242470211330930     ,&
!  0.65234143986191861     ,&
! -0.65599359213247133     ,&
!  -1.5866572642068032E-008,&
!   2.6106713597762773E-002,&
!  0.51345916802183456     ,&
!  -1.1944735457452070E-007,&
!  0.76790299230214776     ,&
!  0.30755434142152160     ,&
! -0.19309092509172471     ,&
!   7.3991071178379259E-008,&
! -0.12170003659699385     ,&
!  0.57319305435412815     ,&
!  0.67858187094366040     ,&
! -0.42437029397046067     ,&
!  0.12137768656752086     ,&
!  -1.3140837692319078E-002,&
!  0.10642233864661235     ,&
!   6.8229626171177837E-002,&
!  0.54172615932574197     ,&
! -0.71799809717192553     ,&
! -0.40107359080949945     ,&
!  0.11471434951865310     ,&
!  -1.2419438850905811E-002,&
! -0.11260400845915516     ,&
!   6.4484016212712075E-002,&
!   6.5930973512194166E-002,&
! -0.11228972412989915     ,&
!   6.8688492559553518E-002,&
! -0.35608684886894670     ,&
! -0.39068637722089328     ,&
!  0.71599413201580542     ,&
!  0.43156080592657275     ,&
!   6.9349962933107787E-002,&
!  0.10675389800335998     ,&
!   7.2250510726294487E-002,&
! -0.37455205443982248     ,&
! -0.41094580114670282     ,&
! -0.68069598649826357     ,&
!  0.45393978962194220     ,&
!  -8.2877165325514465E-002,&
!   1.1385499303904153E-017,&
!  0.11963412299910281     ,&
!  0.42997902589216386     ,&
!  0.45885627469354723     ,&
!   5.9782336257890463E-017,&
!  0.76379840831794654/), (/homo,homo/))     
   print *, "Q: ", Q
   
   ! Convert Q to full matrix type
   nullify(rmat)
   CALL cp_fm_create(rmat,zij(1,1)%matrix%matrix_struct)
   rmat%matrix_struct%nrow_global = homo
   rmat%matrix_struct%ncol_global = homo
   CALL cp_fm_set_all(rmat, 0._dp, 1._dp)
   rmat%local_data = REAL(Q, dp)
   
   ! Calculate Phi
   !Call rotate_orbitals(rmat,vectors)             ! Produces orthogonalized localized orbitals
   vectors%local_data = matmul(real(vectors%local_data),real(Q))
   ! Clean up
   CALL cp_fm_release(matrix=rmat)
   CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
   CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)
   print *, "end scdm calculations"

  END SUBROUTINE scdm_rs

END MODULE scdm
