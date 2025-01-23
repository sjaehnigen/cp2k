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
#include "../base/base_uses.f90"

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
   PUBLIC :: girsanov_type, &
             create_girsanov_type, &
             create_girsanov_environment, &
             release_girsanov_type!, &
!             release_girsanov_environment

! **************************************************************************************************
   TYPE girsanov_type
      INTEGER                          :: n0 = 0, n1 = 0
      ! futher needs: delta_t, Temp, tau, mass_ll of the (selected) particles
      ! INTEGER                          :: virial_components = do_clv_geo_center
      REAL(KIND=dp)                    :: temperature = 0.0_dp
      ! TYPE(npt_info_type), POINTER  :: npt(:, :) => NULL()
      TYPE(section_vals_type), POINTER :: section => NULL()
   END TYPE girsanov_type

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'girsanov_types'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param girsanov ...
!> \param md_section ...
!> \param force_env ...
!> \param simpar ...
!> \param globenv ...
!> \par History
!>      12.2024 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE create_girsanov_type(girsanov, girsanov_section) !, force_env, simpar,  globenv)
      TYPE(girsanov_type), POINTER                       :: girsanov
 !      TYPE(force_env_type), POINTER                      :: force_env
!      TYPE(simpar_type), POINTER                         :: simpar
 !      TYPE(global_environment_type), POINTER             :: globenv

!       LOGICAL                                            :: check, explicit
!  !      TYPE(cell_type), POINTER                           :: cell
      TYPE(section_vals_type), POINTER                   :: girsanov_section
! 
!       check = .NOT. ASSOCIATED(girsanov)
!       CPASSERT(check)
!      girsanov_section => section_vals_get_subs_vals(md_section, "GIRSANOV")
      CALL section_vals_get(girsanov_section) !, explicit=explicit)
!       ! check: ensemble, thermostat, metadynamics or any other bias, selection/colvar
!       ! SEE WHAT HAPPENS WHEN SECTION IS MISSING
!       IF (simpar%ensemble == npt_i_ensemble .OR. &
!           simpar%ensemble == npt_ia_ensemble .OR. &
!           simpar%ensemble == npt_f_ensemble .OR. &
!           simpar%ensemble == npe_f_ensemble .OR. &
!           simpar%ensemble == npe_i_ensemble .OR. &
!           simpar%ensemble == nph_uniaxial_ensemble .OR. &
!           simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
!          ALLOCATE (girsanov)
!          girsanov%section => girsanov_section
!          ! CALL force_env_get(force_env, cell=cell)
! 
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
!       END IF
      NULLIFY (girsanov)
      !girsanov%temperature = 100.0_dp

   END SUBROUTINE create_girsanov_type

! **************************************************************************************************
!> \brief ...
!> \param barostat ...
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   
   SUBROUTINE release_girsanov_type(girsanov)
      TYPE(girsanov_type), POINTER                       :: girsanov
!       TYPE(barostat_type), POINTER                       :: barostat
! 
!       LOGICAL                                            :: check
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
      !girsanov%temperature = 100
      NULLIFY (girsanov)

   END SUBROUTINE release_girsanov_type

! **************************************************************************************************
!> \brief creates the Girsanov environment
!> \param girsanov_env ...
!> \param girsanov_section ...
!> \par History
!>      01.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE create_girsanov_environment(girsanov_env, md_section)
      TYPE(girsanov_type), POINTER                       :: girsanov_env
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

      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: girsanov_section 

      ! CALL cite_reference(VandeVondele2002)

      girsanov_section => section_vals_get_subs_vals(md_section, "GIRSANOV")
      CALL section_vals_get(girsanov_section, explicit=explicit)

      IF (explicit) THEN

          CALL cp_warn(__LOCATION__, &
                         ">>> I was told to create a Girsanov section! "// &
                         "Congratulations!")
          CALL create_girsanov_type(girsanov_env, girsanov_section)

      END IF

   END SUBROUTINE create_girsanov_environment

! ! **************************************************************************************************
! !> \brief releases the fe_env
! !> \param fe_env ...
! !> \author Laino Teodoro 01.2007
! ! **************************************************************************************************
!   SUBROUTINE release_girsanov_environment(girsanov_env, girsanov_section)
!   END SUBROUTINE release_girsanov_environment
!   SUBROUTINE release_thermostat_type(thermostat)
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
!       IF (ASSOCIATED(fe_env)) THEN
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
!    END SUBROUTINE fe_env_release

END MODULE girsanov_types
