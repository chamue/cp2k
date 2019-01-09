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
   USE qs_environment_types,             ONLY: get_qs_env, &
                                               qs_environment_type
   USE atomic_kind_types,                ONLY: atomic_kind_type
   USE qs_kind_types,                    ONLY: qs_kind_type
   USE cell_types,                       ONLY: cell_type
   USE cp_control_types,                 ONLY: dft_control_type
   USE qs_mo_types,                      ONLY: mo_set_p_type
   USE particle_types,                   ONLY: particle_type
   USE pw_env_types,                     ONLY: pw_env_get,&
                                               pw_env_type    
   USE pw_pool_types,                    ONLY: pw_pool_type, &
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
   USE lapack,                          ONLY: lapack_dgeqp3
   USE kinds,                           ONLY: dp
   USE pw_grid_types,                   ONLY: PW_MODE_LOCAL
   USE message_passing,                 ONLY: mp_recv, &
                                              mp_send, &
                                              mp_sum, &
                                              mp_sync

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
  SUBROUTINE scdm_rs(qs_env, zij, vectors, ispin)
   ! Input Variables
   TYPE(qs_environment_type), POINTER            :: qs_env
   TYPE(cp_fm_p_type), INTENT(INOUT)             :: zij(:,:)
   TYPE(cp_fm_type), POINTER                     :: vectors
   INTEGER                                       :: ispin

   ! Working Variables
   integer                                       :: k, i, j, INFO, reason, M, N, &
                                                    ivector, homo
   real(kind=dp), allocatable, dimension(:,:)    :: vectors_rs_2d, Q, H, Vtensor
   real(kind=dp), allocatable, dimension(:)      :: Tau, v, Work
   real(kind=dp)                                 :: omega
   integer, allocatable, dimension(:)            :: jpvt 
   TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                    :: atomic_kind_set !!! get_qs_env
   TYPE(cp_fm_type), POINTER                     :: mo_coeff        !!! get_mo_set
   TYPE(dft_control_type), POINTER               :: dft_control     !!! get_qs_env
   TYPE(mo_set_p_type), DIMENSION(:), &
       POINTER                                   :: mos             !!! get_qs_env
   TYPE(particle_type), DIMENSION(:), &
      POINTER                                    :: particle_set    !!! get_qs_env
   TYPE(pw_env_type), POINTER                    :: pw_env          !!! get_qs_env
   TYPE(pw_p_type)                               :: wf_g, wf_r      
   TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                    :: qs_kind_set     !!! get_qs_env
   TYPE(pw_pool_type), POINTER                   :: auxbas_pw_pool  !!! pw_env_get
   TYPE(cell_type), POINTER                      :: cell            !!! get_qs_env
   LOGICAL                                       :: failure
   type(cp_fm_type), pointer                     :: rmat           
   integer, dimension(3)                         :: L, U
   CHARACTER(len=*), PARAMETER :: routineN = 'scdm_rs', &
           routineP = moduleN//':'//routineN

   ! Get working variables
   NULLIFY(dft_control, pw_env, auxbas_pw_pool, mos, mo_coeff)
   CALL get_qs_env(qs_env=qs_env,dft_control=dft_control,mos=mos, pw_env=pw_env)
   CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
   CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE)
   ! When no wf_g, error 'Missing actual argument' in calculate_wavefunction
   CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE)
   CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
   failure=.FALSE.

   ! Looping over occupied orbitals and inserting real-space values in matrix columns
   DO ivector=1,homo
         CALL get_qs_env(qs_env=qs_env,&
                         atomic_kind_set=atomic_kind_set,&
                         qs_kind_set=qs_kind_set,&
                         cell=cell,&
                         particle_set=particle_set)
         CALL calculate_wavefunction(mo_vectors=mo_coeff,ivector=ivector,rho=wf_r,&
                 rho_gspace=wf_g,atomic_kind_set=atomic_kind_set,qs_kind_set=qs_kind_set,&
              cell=cell,dft_control=dft_control,particle_set=particle_set,&
              pw_env=pw_env)
          !omega used for normalizing the matrix later on
          omega = wf_r%pw%pw_grid%dh(1,1)*2
          omega = omega*omega*omega
          omega = sqrt(omega)
          if (ivector == 1) then
             !only works for serial calculation, otherwise other upper and lower boundaries
             L(1) = wf_r%pw%pw_grid%bounds(1, 1)
             L(2) = wf_r%pw%pw_grid%bounds(1, 2)
             L(3) = wf_r%pw%pw_grid%bounds(1, 3)
             U(1) = wf_r%pw%pw_grid%bounds(2, 1)
             U(2) = wf_r%pw%pw_grid%bounds(2, 2)
             U(3) = wf_r%pw%pw_grid%bounds(2, 3)
            ALLOCATE(vectors_rs_2d(homo,((U(1)-L(1))/2+1)*((U(2)-L(2))/2+1)*((U(3)-L(1))/2+1)))
          END IF
          call pw_to_matrix(wf_r%pw, vectors_rs_2d, ivector)
   end do

   ! Normalize orbitals
   vectors_rs_2d = vectors_rs_2d*omega

   ! Allocating work variables
   N = size(vectors_rs_2d,dim=1)
   M = size(vectors_rs_2d,dim=2)
   allocate(Tau(N))
   allocate(v(N))
   allocate(Work(3*M+1))
   allocate(Q(N,N))
   allocate(H(N,N))
   allocate(jpvt(M))
   allocate(Vtensor(N,N))
   ! Pivoted QR Decomposition of Psicon
   jpvt = 0
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
   
   ! Convert Q to full matrix type
   nullify(rmat)
   CALL cp_fm_create(rmat,zij(1,1)%matrix%matrix_struct)
   rmat%matrix_struct%nrow_global = homo
   rmat%matrix_struct%ncol_global = homo
   CALL cp_fm_set_all(rmat, 0._dp, 1._dp)
   rmat%local_data = REAL(Q, dp)
   
   ! Calculate Phi
   vectors%local_data = matmul(real(vectors%local_data),real(Q))
   ! Clean up
   CALL cp_fm_release(matrix=rmat)
   CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
   CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)

  END SUBROUTINE scdm_rs

SUBROUTINE pw_to_matrix(pw, matrix, ivector)
      ! Nearly identical to pw_to_cube except for the printing to a matrix instead
      ! of a cube file (and therefore obviously no header)
      TYPE(pw_type), POINTER                             :: pw
      real(kind=dp), dimension(:,:)                      :: matrix
      integer, intent(in)                                :: ivector

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_to_matrix', routineP = moduleN//':'//routineN

      INTEGER :: checksum, dest, gid, handle, I1, I2, I3, ip, L1, L2, L3, my_rank, &
         my_stride(3), num_pe, rank(2), source, tag, U1, U2, U3, c
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: buf
      CALL timeset(routineN, handle)
      my_stride = 2
 
      ! shortcut
      L1 = pw%pw_grid%bounds(1, 1)
      L2 = pw%pw_grid%bounds(1, 2)
      L3 = pw%pw_grid%bounds(1, 3)
      U1 = pw%pw_grid%bounds(2, 1)
      U2 = pw%pw_grid%bounds(2, 2)
      U3 = pw%pw_grid%bounds(2, 3)

      ALLOCATE (buf(L3:U3))

      my_rank = pw%pw_grid%para%my_pos
      gid = pw%pw_grid%para%group
      num_pe = pw%pw_grid%para%group_size
      tag = 1
      rank(2) = my_rank
      checksum = 1
      CALL mp_sum(checksum, gid)
      CPASSERT(checksum == 1)
      dest = rank(2)
      c=1
      DO I1 = L1, U1, my_stride(1)
         DO I2 = L2, U2, my_stride(2)

               ! cycling through the CPUs, check if the current ray (I1,I2) is local to that CPU
               IF (pw%pw_grid%para%mode .NE. PW_MODE_LOCAL) THEN
                  DO ip = 0, num_pe-1
                     IF (pw%pw_grid%para%bo(1, 1, ip, 1) <= I1-L1+1 .AND. pw%pw_grid%para%bo(2, 1, ip, 1) >= I1-L1+1 .AND. &
                         pw%pw_grid%para%bo(1, 2, ip, 1) <= I2-L2+1 .AND. pw%pw_grid%para%bo(2, 2, ip, 1) >= I2-L2+1) THEN
                        source = ip
                     ENDIF
                  ENDDO
               ELSE
                  source = dest
               ENDIF

               IF (source == dest) THEN
                  IF (my_rank == source) THEN
                     buf(:) = pw%cr3d(I1, I2, :)
                  ENDIF
               ELSE
                  IF (my_rank == source) THEN
                     buf(:) = pw%cr3d(I1, I2, :)
                     CALL mp_send(buf, dest, tag, gid)
                  ENDIF
                  IF (my_rank == dest) THEN
                     CALL mp_recv(buf, source, tag, gid)
                  ENDIF
               ENDIF
               do I3 = L3,U3,my_stride(3)
                  matrix(ivector,c) = buf(I3)
                  c = c+1 
               end do

               ! this double loop generates so many messages that it can overload
               ! the message passing system, e.g. on XT3
               ! we therefore put a barrier here that limits the amount of message
               ! that flies around at any given time.
               ! if ever this routine becomes a bottleneck, we should go for a
               ! more complicated rewrite
               CALL mp_sync(gid)

         ENDDO
      ENDDO

      DEALLOCATE (buf)

      CALL timestop(handle)

   END SUBROUTINE pw_to_matrix
end module scdm

