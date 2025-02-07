!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types and their creation routines for Girsanov reweighting of an MD trajectory
!> \par History
!>      12.2024 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************

MODULE girsanov_types
!   USE global_types,                    ONLY: global_environment_type
   USE input_constants,                 ONLY: grw_bias,&
                                              grw_perturbation
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type!,&
!                                               section_vals_val_get
   USE kinds,                           ONLY: dp
   USE simpar_types,                    ONLY: simpar_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

!   INTEGER, PARAMETER, PUBLIC               :: do_clv_geo_center = 0, &
!                                               do_clv_fix_point = 1, &
!                                               do_clv_xyz = 0, &
!                                               do_clv_x = 1, &
!                                               do_clv_y = 2, &
!                                               do_clv_z = 3, &
!                                               do_clv_xy = 4, &
!                                               do_clv_xz = 5, &
!                                               do_clv_yz = 6

   PRIVATE
   PUBLIC :: girsanov_env_type, &
             rw_factor_type, &
             create_girsanov_environment, &
             release_girsanov_environment

! **************************************************************************************************
!> \brief defines meta_env type
! **************************************************************************************************
!   TYPE meta_env_type
!      LOGICAL                                   :: do_hills = .FALSE., do_multiple_walkers = .FALSE.
!      LOGICAL                                   :: extended_lagrange = .FALSE.
!      LOGICAL                                   :: well_tempered = .FALSE.
!      LOGICAL                                   :: langevin = .FALSE.
!      LOGICAL                                   :: use_plumed = .FALSE.
!      CHARACTER(LEN=default_path_length)        :: plumed_input_file = ""
!      INTEGER                                   :: n_colvar = -1
!      REAL(KIND=dp)                             :: ekin_s = 0.0_dp, epot_s = 0.0_dp, dt = 0.0_dp, &
!                                                   avg_temp = 0.0_dp, epot_walls = 0.0_dp
!      LOGICAL                                   :: tempcontrol = .FALSE., restart = .FALSE.
!      REAL(KIND=dp)                             :: temp_wanted = 0.0_dp, toll_temp = 0.0_dp
!      REAL(KIND=dp)                             :: delta_t = 0.0_dp, invdt = 0.0_dp, &
!                                                   wtgamma = 0.0_dp, wttemperature = 0.0_dp
!      INTEGER                                   :: n_steps = -1
!      ! time pointer should never be allocated itself.. that's a quite bad
!      ! choice derived from the md_env.. So whenever the meta_env type is
!      ! used the pointer time should be pointer to the high level time counter
!      REAL(KIND=dp), POINTER                    :: time => NULL()
!      TYPE(hills_env_type), POINTER             :: hills_env => NULL()
!      TYPE(metavar_type), POINTER, DIMENSION(:) :: metavar => NULL()
!      TYPE(multiple_walkers_type), POINTER      :: multiple_walkers => NULL()
!      TYPE(mp_para_env_type), POINTER           :: para_env => NULL()
!      TYPE(section_vals_type), POINTER          :: metadyn_section => NULL()
!      TYPE(rng_stream_type), DIMENSION(:), &
!         ALLOCATABLE :: rng
!      INTEGER                                   :: TAMCSteps = -1
!      REAL(kind=dp)                             :: zdt = 0.0_dp
!   END TYPE meta_env_type



! **************************************************************************************************
!> \brief defines the environment type for Girsanov reweighting
! **************************************************************************************************
   TYPE girsanov_env_type
      REAL(KIND=dp)                                      :: n0 = 0, n1 = 0
      REAL(KIND=dp)                                      :: fx = 0, fy = 0, fz=0
      TYPE(rw_factor_type), POINTER, DIMENSION(:)        :: rw_factor => NULL()
      TYPE(section_vals_type), POINTER                   :: girsanov_section => NULL()

      ! futher needs: delta_t, Temp, tau, mass_ll of the (selected) particles
      ! INTEGER                                        :: virial_components = do_clv_geo_center
      ! REAL(KIND=dp)                                  :: temperature = 0.0_dp
!       CHARACTER(LEN=default_path_length), &
!          DIMENSION(:), POINTER                :: walkers_file_name => NULL()

   END TYPE girsanov_env_type

! **************************************************************************************************
!> \brief defines the type for the reweighting factor and buffering
! **************************************************************************************************
   TYPE rw_factor_type
      INTEGER                                            :: n0 = 0, n1 = 0
      ! futher needs: delta_t, Temp, tau, mass_ll of the (selected) particles
      ! INTEGER                          :: virial_components = do_clv_geo_center
   END TYPE rw_factor_type

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'girsanov_types'

CONTAINS


! **************************************************************************************************
!> \brief creates the Girsanov environment
!> \param girsanov_env ...
!> \param girsanov_section ...
!> \par History
!>      01.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE create_girsanov_environment(girsanov_env, md_section)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      ! TYPE(section_vals_type), POINTER                   :: fe_section
   
      ! INTEGER                                            :: i, id_method
      ! LOGICAL                                            :: explicit
      ! TYPE(section_vals_type), POINTER                   :: ui_section, ui_var_section

      ! CPASSERT(.NOT. ASSOCIATED(fe_env))
      !    
      ! CALL section_vals_get(fe_section, explicit=explicit)
      ! IF (explicit) THEN
      ! END IF
      TYPE(section_vals_type), POINTER                   :: md_section

      LOGICAL                                            :: check, explicit
      TYPE(section_vals_type), POINTER                   :: girsanov_section 

      ! CALL cite_reference(VandeVondele2002)

      check = .NOT. ASSOCIATED(girsanov_env)
      CPASSERT(check)

      girsanov_section => section_vals_get_subs_vals(md_section, "GIRSANOV")
      CALL section_vals_get(girsanov_section, explicit=explicit)

      IF (explicit) THEN

          CALL cp_warn(__LOCATION__, &
                         ">>> I was told to create a Girsanov section! "// &
                         "Congratulations!")
          ! CALL create_girsanov_type(girsanov_env, girsanov_section)
          ALLOCATE (girsanov_env)
          girsanov_env%n0 = 100_dp

      END IF
!          girsanov%temperature = simpar%temp_ext
!          ! CALL section_vals_val_get(girsanov_section, "TEMP_TOL", r_val=simpar%temp_baro_tol)
!          ! ! Initialize or possibly restart Barostat
!          ! CALL initialize_npt(simpar, globenv, barostat%npt, &
!          !                     cell, work_section=barostat_section)
! 
!          ! User defined virial screening
!          CALL section_vals_val_get(girsanov_section, "VIRIAL", i_val=barostat%virial_components)
!          check = barostat%virial_components == do_clv_xyz .OR. simpar%ensemble == npt_f_ensemble
!          IF (.NOT. check) &
!             CALL cp_abort(__LOCATION__, "The screening of the components of "// &
!                           "the virial is available only with the NPT_F ensemble!")
!       ELSE
!          IF (explicit) &
!             CALL cp_warn(__LOCATION__, &
!                          "A barostat has been defined with an MD ensemble which does not support barostat! "// &
!                          "Its definition will be ignored!")

   END SUBROUTINE create_girsanov_environment

! **************************************************************************************************
!> \brief releases the Girsanov environment
!> \param girsanov_env ...
!> \par History
!>      01.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE release_girsanov_environment(girsanov_env)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env

!      TYPE(thermostat_type), INTENT(INOUT)               :: thermostat
!
!      NULLIFY (thermostat%section)
!      IF (ASSOCIATED(thermostat%nhc)) THEN
!         CALL lnhc_dealloc(thermostat%nhc)
!      END IF
!      IF (ASSOCIATED(thermostat%csvr)) THEN
!         CALL csvr_dealloc(thermostat%csvr)
!      END IF
! 
!      OR

!    SUBROUTINE fe_env_release(fe_env)
!       TYPE(free_energy_type), POINTER                    :: fe_env
! 
!       INTEGER                                            :: i
! 
! !      IF (ASSOCIATED(barostat)) THEN
!          check = barostat%ref_count > 0
!          CPASSERT(check)
!          barostat%ref_count = barostat%ref_count - 1
!          IF (barostat%ref_count < 1) THEN
!             IF (ASSOCIATED(barostat%npt)) THEN
!                DEALLOCATE (barostat%npt)
!             END IF
!             NULLIFY (barostat%section)
!             DEALLOCATE (barostat)
!          END IF
!       END IF
      IF (ASSOCIATED(girsanov_env)) THEN
          ! NULLIFY SECTION
          ! CALL release_girsanov_type(girsanov_env)
          DEALLOCATE (girsanov_env)
      END IF
!       IF (ASSOCIATEDV(fe_env)) THEN
!          IF (ASSOCIATED(fe_env%covmx)) THEN
!             DEALLOCATE (fe_env%covmx)
!          END IF
!          IF (ASSOCIATED(fe_env%cg_data)) THEN
!             DO i = 1, SIZE(fe_env%cg_data)
!                IF (ASSOCIATED(fe_env%cg_data(i)%avg)) THEN
!                   DEALLOCATE (fe_env%cg_data(i)%avg)
!                END IF
!                IF (ASSOCIATED(fe_env%cg_data(i)%var)) THEN
!                   DEALLOCATE (fe_env%cg_data(i)%var)
!                END IF
!             END DO
!             DEALLOCATE (fe_env%cg_data)
!          END IF
!          IF (ASSOCIATED(fe_env%conv_par)) THEN
!             DEALLOCATE (fe_env%conv_par)
!          END IF
!          IF (ASSOCIATED(fe_env%uivar)) THEN
!             DO i = 1, SIZE(fe_env%uivar)
!                IF (ASSOCIATED(fe_env%uivar(i)%ss)) THEN
!                   DEALLOCATE (fe_env%uivar(i)%ss)
!                END IF
!             END DO
!             DEALLOCATE (fe_env%uivar)
!          END IF
!          DEALLOCATE (fe_env)
!       END IF
   END SUBROUTINE release_girsanov_environment

END MODULE girsanov_types
