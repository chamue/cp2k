!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Localization method scdm using pivoted QR decomposition
!> \par History Implementation scdm during master thesis
!>
!> \author Charlotte Mueller (11.2018) 
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
   USE particle_list_types,              ONLY: particle_list_create, &
                                               particle_list_type
   USE pw_env_types,                     ONLY: pw_env_get,&
                                               pw_env_type    
   USE pw_pool_types,                    ONLY: pw_pool_type, &
                                               pw_pool_p_type, &
                                               pw_pool_create_pw
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
   USE message_passing,                  ONLY: mp_maxloc,&
                                               mp_recv,&
                                               mp_send,&
                                               mp_sum, &
                                               mp_sync
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
   integer                                       :: k, l, i, j, INFO, reason, M, N, &
                                                    ivector, ne, np, ii1, ii2, ii3
   real(kind=dp), allocatable, dimension(:,:)    :: vectors_rs ! M, N
   double precision, allocatable, dimension(:,:) :: Psicon ! N, M
   double precision, allocatable, dimension(:)   :: Tau, v ! N
   double precision, allocatable, dimension(:)   :: Work ! 3M+1
   double precision, allocatable, dimension(:,:) :: Q, H, Vtensor ! N, N
   double precision                              :: omega
   integer, allocatable, dimension(:)            :: jpvt 
   INTEGER :: homo, my_stride
   TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                    :: atomic_kind_set
   TYPE(cp_fm_type), POINTER                     :: mo_coeff
   TYPE(dft_control_type), POINTER               :: dft_control
   TYPE(mo_set_p_type), DIMENSION(:), &
       POINTER                                   :: mos
   TYPE(particle_list_type), POINTER             :: particles 
   TYPE(particle_type), DIMENSION(:), &
      POINTER                                    :: particle_set
   TYPE(pw_env_type), POINTER                    :: pw_env
   TYPE(pw_p_type)                               :: wf_g, wf_r
   TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                    :: qs_kind_set
   TYPE(particle_list_type), POINTER             :: my_particles  
   TYPE(pw_pool_type), POINTER                   :: auxbas_pw_pool
   TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                    :: pw_pools
   TYPE(cell_type), POINTER                      :: cell
   INTEGER, allocatable, DIMENSION(:)            :: particles_z
   INTEGER :: checksum, dest, gid, I1, I2, I3, N2, N3, ind, ip, L1, L2, &
      L3, my_rank, num_pe, rank(2), source, stat, tag, U1, &
     U2, U3
   LOGICAL                                       :: failure
   REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)      :: buf
   type(cp_fm_type), pointer                     :: rmat

   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN
   
   print *, "I arrived in scdm module"
   ! Start calculating realspace orbitals
   NULLIFY(atomic_kind_set, qs_kind_set, cell, dft_control, mos, &
           particle_set, pw_env)

   CALL get_qs_env(qs_env=qs_env,&
                         atomic_kind_set=atomic_kind_set,&
                         qs_kind_set=qs_kind_set,&
                         cell=cell,&
                         dft_control=dft_control,&
                         mos=mos,&
                         particle_set=particle_set,&
                         pw_env=pw_env)
   print *, "Got the qs environment"
   Nullify(auxbas_pw_pool, pw_pools)
   CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                          pw_pools=pw_pools)
   CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE)
   ! When no wf_g, error 'Missing actual argument' in calculate_wavefunction
   CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE)
   print *, "Created the pw environment"
   CALL particle_list_create(list=particles,els_ptr=particle_set)
   Nullify(mo_coeff)
   CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
   NULLIFY(my_particles)
   my_particles=>particles
   ne = my_particles%n_els
   ALLOCATE(particles_z(ne))
   DO i=1, ne
      CALL get_atomic_kind(my_particles%els(i)%atomic_kind,z=particles_z(i)) 
   END DO
   print *, "Got all the MO information"
   failure=.FALSE.
   my_stride = 1
   ! for now, only stride = 1 is allowed
   !if (present(stride)) then 
   !   CALL cp_assert(SIZE(stride)==1.OR.SIZE(stride)==3,cp_fatal_level,cp_assertion_failed,routineP,&
   !               "STRIDE keyword can accept only 1 (the same for X,Y,Z) or 3 values. Correct your input file."//&
   !               CPSourceFileRef,&
   !               only_ionode=.TRUE.)
   !   IF (SIZE(stride)==3) THEN
   !      DO i = 1,3
   !        CALL cp_assert(stride(1)==stride(2).AND.stride(1)==stride(3),cp_fatal_level,cp_assertion_failed,routineP,&
   !                    "scdm method only available with 3 equal stride values."//&
   !                    CPSourceFileRef,&
   !                    only_ionode=.TRUE.)
   !        my_stride = stride(1)
   !      END DO
   !   ELSE 
   !           my_stride = stride
   !   END IF
   !end if
   !CPPrecondition(my_stride>0,cp_failure_level,routineP,error,failure)
   np=MIN(99999,SIZE(particles_z))
   print *, "Processed Stride"
   ! getting the pw grid
   DO ivector=1,homo
         print *, "at orbital", ivector
         CALL calculate_wavefunction(mo_vectors=mo_coeff,ivector=ivector,rho=wf_r,&
                 rho_gspace=wf_g,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
              cell=cell,dft_control=dft_control,particle_set=particle_set,&
              pw_env=pw_env)
          print *, "calculated wavefunction"
          omega = wf_r%pw%pw_grid%dh(1,1)*REAL(my_stride,dp)
          omega = omega*omega*omega
          omega = sqrt(omega)
          print *, "calculated omega"
          ! shortcut
          L1=wf_r%pw%pw_grid%bounds(1,1)
          L2=wf_r%pw%pw_grid%bounds(1,2)
          L3=wf_r%pw%pw_grid%bounds(1,3)
          U1=wf_r%pw%pw_grid%bounds(2,1)
          U2=wf_r%pw%pw_grid%bounds(2,2)
          U3=wf_r%pw%pw_grid%bounds(2,3)
          ALLOCATE(buf(L3:U3),stat=stat)
          IF (ivector==1) then
                  ALLOCATE(vectors_rs((U1-L1)*(U2-L2)*(U3-L3),homo))
          END IF
          print *, "about to initiate parallel stuff"
          my_rank     = wf_r%pw%pw_grid%para%my_pos
          gid         = wf_r%pw%pw_grid%para%group
          num_pe      = wf_r%pw%pw_grid%para%group_size
          tag         = 1
          print *, "done"
          rank(2)=my_rank
          checksum=1
          print *, "about to call mp functions"
          CALL mp_sum(checksum,gid)
          print *, "checksum: ", checksum
          CPASSERT(checksum==1)
          CALL mp_maxloc(rank,gid)

          print *, "about to start double loop"
          dest = rank(2)
          ii1 = 0
          ii2 = 0
          ii3 = 0
          DO I1=L1,U1,my_stride
            ii1 = ii1 + 1
            DO I2=L2,U2,my_stride
                ii2 = ii2 + 1
                ! cycling through the CPUs, check if the current ray (I1,I2) is local to that CPU
                IF (wf_r%pw%pw_grid%para%mode .NE. PW_MODE_LOCAL ) THEN
                   DO ip=0,num_pe-1
                      IF (wf_r%pw%pw_grid%para%bo (1,1,ip,1)<=I1-L1+1 .AND. wf_r%pw%pw_grid%para%bo (2,1,ip,1)>=I1-L1+1 .AND. &
                          wf_r%pw%pw_grid%para%bo (1,2,ip,1)<=I2-L2+1 .AND. wf_r%pw%pw_grid%para%bo (2,2,ip,1)>=I2-L2+1 ) THEN
                         source = ip
                      ENDIF
                    ENDDO
                ELSE
                   source = dest
                ENDIF
                IF (source==dest) THEN
                   IF (my_rank==source) THEN
                      buf(:)=wf_r%pw%cr3d(I1,I2,:)
                   ENDIF
                ELSE
                   IF(my_rank==source) THEN
                      buf(:)=wf_r%pw%cr3d(I1,I2,:)
                      CALL mp_send(buf,dest,tag,gid)
                   ENDIF
                   IF(my_rank==dest) THEN
                   CALL mp_recv(buf,source,tag,gid)
                   ENDIF
                ENDIF
                DO I3 = L3,U3,my_stride
                   ii3 = ii3 + 1
                   N2 = U2 - L2
                   N3 = U3 - L3
                   ind = ((ii1-1)*N2+ii2-1)*N3+ii3
                   vectors_rs(ind,ivector)=buf(I3)
                END DO

          ! this double loop generates so many messages that it can overload
          ! the message passing system, e.g. on XT3
          ! we therefore put a barrier here that limits the amount of message
          ! that flies around at any given time.
          ! if ever this routine becomes a bottleneck, we should go for a
          ! more complicated rewrite
               CALL mp_sync(gid)

         ENDDO
       ENDDO
       
       DEALLOCATE(buf)
   END DO
   
   print *, "start scdm calculation"
   M = size(vectors_rs,dim=1)
   N = size(vectors_rs,dim=2)
   ! Allocate Arrays
   allocate(Psicon(N, M))
   allocate(Tau(N))
   allocate(v(N))
   allocate(Work(3*M+1))
   allocate(Q(N,N))
   allocate(H(N,N))
   allocate(jpvt(M))
   allocate(Vtensor(N,N))

   vectors_rs = vectors_rs*omega

   ! Calculate Psicon
   Psicon = TRANSPOSE(DBLE(vectors_rs))

   ! Pivoted QR Decomposition of Psicon
   call lapack_dgeqp3(N,M,Psicon,N,jpvt,Tau,Work,3*M+1,INFO)
   print *, 'INFO (0 = success)', INFO
do i = 1,N
      H = 0
      v(1:i-1)=0
      v(i)=1
      v(i+1:N)=Psicon(i+1:N,i)
      do j = 1,N
         H(j,j)=1
      end do
      do k = 1,N
        do l = 1,N
          Vtensor(k,l) = v(k)*v(l)
        end do
      end do
      H = H - Tau(i)*Vtensor
      if (i == 1) then
         Q = H
      else
         Q = matmul(Q,H)
      end if
   end do
   print *, "Created Q (dp) successfully"
   ! Convert Q to full matrix type
   nullify(rmat)
   CALL cp_fm_create(rmat,zij(1,1)%matrix%matrix_struct)
   rmat%matrix_struct%nrow_global = homo
   rmat%matrix_struct%ncol_global = homo
   print *, "Initialized rmat successfully"
   CALL cp_fm_set_all(rmat, 0._dp, 1._dp)
   print *, "Set rmat"
   rmat%local_data = REAL(Q, dp)
   print *, "Inserted Q in rmat"
   ! Calculate Phi
   Call rotate_orbitals(rmat,vectors)                                ! Produces orthogonalized Phi
   print *, "Rotated orbitals"
   !vectors = matmul(vectors,transpose(vectors(jpvt(1:N),:))) ! non-orthogonal Phi
   CALL cp_fm_release(matrix=rmat)
   CALL cp_fm_release(matrix=mo_coeff)
   print *, "released rmat"
   deallocate(Psicon, Tau, v, Work, Q, H, jpvt, Vtensor, vectors_rs, particles_z) 
   print *, "end"

  END SUBROUTINE scdm_rs

END MODULE scdm

