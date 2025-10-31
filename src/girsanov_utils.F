!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Girsanov utils
!> \par History
!>      12.2024 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
MODULE girsanov_utils
   USE girsanov_types,                  ONLY: girsanov_env_type,&
                                              roll_girsanov_environment,&
                                              degree_of_freedom_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: grw_bias,&
                                              grw_perturbation
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              default_string_length
   USE machine,                         ONLY: m_flush
   USE message_passing,                 ONLY: mp_comm_type   
   USE physcon,                         ONLY: massunit
                                              
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: calculate_reweighting_factor, &
             write_girsanov_output
             

! *** Global parameters ***
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'girsanov_utils'

CONTAINS


! **************************************************************************************************
!> \brief Calculates the Girsanov reweighting factor(s) based on integration step data
!> \param girsanov_env ...
!> \param group ...
!> \par History
!>      06.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE calculate_reweighting_factor(girsanov_env, group)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      CLASS(mp_comm_type), INTENT(IN)                    :: group

      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: dof
      REAL(KIND=dp), POINTER                             :: S_k
      REAL(KIND=dp)                                      ::  weight
      INTEGER                                            :: j

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'calculate_reweighting_factor'
      INTEGER                                            :: handle
      
      CALL timeset(routineN, handle)

      dof => girsanov_env%degrees_of_freedom
      S_k => girsanov_env%delta_S

      IF (ASSOCIATED(girsanov_env)) THEN

         CALL group%sum(dof(:)%eta(1))
         CALL group%sum(dof(:)%eta(2))
         ! NB: para_sum for grad_U done by force_env!

         DO j=1, SIZE(dof)
            weight = SQRT(dof(j)%mass)
            dof(j)%delta_eta(1) = girsanov_env%prefactor_eta_1 / weight * dof(j)%grad_U(1)
            dof(j)%delta_eta(2) = girsanov_env%prefactor_eta_2 / weight * dof(j)%grad_U(2)
            S_k = S_k + dof(j)%eta(1)*dof(j)%delta_eta(1) + 0.5_dp*dof(j)%delta_eta(1)**2 &
                      + dof(j)%eta(2)*dof(j)%delta_eta(2) + 0.5_dp*dof(j)%delta_eta(2)**2
         ENDDO
      END IF

      NULLIFY(dof)

      CALL timestop(handle)

   END SUBROUTINE calculate_reweighting_factor


! **************************************************************************************************
!> \brief Parser for Girsanov output
!> \param girsanov_env ...
!> \param motion_section ...
!> \param itimes ...
!> \param time ...
!> \param pos ...
!> \param act ....
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE write_girsanov_output(girsanov_env, motion_section, itimes, time, pos, act)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(section_vals_type), POINTER                   :: motion_section
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      CHARACTER(LEN=default_string_length), INTENT(IN), &
         OPTIONAL                                        :: pos, act

      CHARACTER(LEN=default_string_length)               :: my_act, my_pos

      IF (ASSOCIATED(girsanov_env)) THEN
         my_pos = "APPEND"
         my_act = "WRITE"
         IF (PRESENT(pos)) my_pos = pos
         IF (PRESENT(act)) my_act = act

         CALL write_girsanov_output_low(girsanov_env, motion_section, itimes, time, my_pos, my_act, &
                                        "BIAS_FORCE", ".girsanov_bias")
         CALL write_girsanov_output_low(girsanov_env, motion_section, itimes, time, my_pos, my_act, &
                                        "RANDOM_NUMBERS", ".girsanov_eta")
         CALL write_girsanov_output_low(girsanov_env, motion_section, itimes, time, my_pos, my_act, &
                                        "REWEIGHTING_FACTOR", ".girsanov_factor")
      END IF

   END SUBROUTINE write_girsanov_output

! **************************************************************************************************
!> \brief Parser for Girsanov output
!> \param girsanov_env ...
!> \param motion_section ...
!> \param itimes ...
!> \param time ...
!> \param pos ...
!> \param act ....
!> \param key_suffix ....
!> \param extension ....
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE write_girsanov_output_low(girsanov_env, motion_section, itimes, time, pos, act, &
                                        key_suffix, extension)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(section_vals_type), POINTER                   :: motion_section
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      CHARACTER(*), INTENT(IN)                           :: pos, act, key_suffix, extension

      INTEGER                                            :: output_unit
      CHARACTER(LEN=default_string_length)               :: full_key
      LOGICAL                                            :: new_file
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      full_key = "MD%GIRSANOV%PRINT%" // TRIM(key_suffix)
      ! Print Bias Force (DEBUG)
      output_unit = cp_print_key_unit_nr(logger, motion_section, full_key, &
                                         extension=TRIM(extension), file_position=pos, &
                                         file_action=act, file_form="FORMATTED", &
                                         is_new_file=new_file)
      IF (output_unit > 0) THEN    
         SELECT CASE (key_suffix)
         CASE ("BIAS_FORCE")
            CALL write_bias_force(girsanov_env, itimes, time, output_unit)
         CASE ("RANDOM_NUMBERS")
            CALL write_random_numbers(girsanov_env, itimes, time, output_unit)
         CASE ("REWEIGHTING_FACTOR")
            CALL write_reweighting_factor(girsanov_env, itimes, time, output_unit, new_file)
         END SELECT
      END IF

      CALL cp_print_key_finished_output(output_unit, logger, motion_section, full_key)

   END SUBROUTINE write_girsanov_output_low

! **************************************************************************************************
!> \brief Write Girsanov bias force per DOF during an MD run
!> \param girsanov_env ...
!> \param itimes ...
!> \param time ...
!> \param output_unit ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE write_bias_force(girsanov_env, itimes, time, output_unit)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: dof      
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: i, n, n_atoms, n_dim, i_entry, i_particle
      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map
      INTEGER, POINTER, DIMENSION(:)                     :: atom_list

      dof => girsanov_env%degrees_of_freedom
      n_atoms = girsanov_env%n_atoms
      n_dim = girsanov_env%n_dim
      atom_list => girsanov_env%atom_list
      index_map => girsanov_env%index_map

      WRITE (output_unit, "(I8)") n_atoms
      ! WRITE (output_unit, "(A,I8,A,F12.3,A,F20.10)") &
      WRITE (output_unit, "(A,I8,A,F12.3,A,ES25.16E2)") &
         " i = ", itimes, ", time = ", time, ", U = ", girsanov_env%U
      DO i = 1, n_atoms
         DO  n = 1, n_dim
            i_particle = atom_list(i)
            i_entry = index_map(n, i_particle) 
            WRITE (output_unit, "(T2,I8,I8,1X,F15.6,1X,ES25.16E2)") &
               i_particle, n, dof(i_entry)%mass/massunit, dof(i_entry)%grad_U(2)
         END DO
      END DO

      CALL m_flush(output_unit)
            
      NULLIFY(index_map)
      NULLIFY(atom_list)

   END SUBROUTINE write_bias_force

! **************************************************************************************************
!> \brief Write Girsanov random numbers per DOF during an MD run
!> \param girsanov_env ...
!> \param itimes ...
!> \param time ...
!> \param output_unit ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE write_random_numbers(girsanov_env, itimes, time, output_unit)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: dof      
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: i, n, n_atoms, n_dim, i_entry, i_particle
      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map
      INTEGER, POINTER, DIMENSION(:)                     :: atom_list

      dof => girsanov_env%degrees_of_freedom
      n_atoms = girsanov_env%n_atoms
      n_dim = girsanov_env%n_dim
      atom_list => girsanov_env%atom_list
      index_map => girsanov_env%index_map

      WRITE (output_unit, "(I8)") n_atoms
      WRITE (output_unit, "(A,I8,A,F12.3,A,F20.10)") &
         " i = ", itimes, ", time = ", time, ", U = ", girsanov_env%U
      DO i = 1, n_atoms
         DO  n = 1, n_dim
            i_particle = atom_list(i)
            i_entry = index_map(n, i_particle)             
            WRITE (output_unit, "(T2,I8,I8,4(1X,ES25.16E2))") &
               i_particle, n, dof(i_entry)%eta(1:2), dof(i_entry)%delta_eta(1:2)
         END DO
      END DO

      CALL m_flush(output_unit)

      NULLIFY(index_map)
      NULLIFY(atom_list)   

   END SUBROUTINE write_random_numbers

! **************************************************************************************************
!> \brief Write Girsanov reweighting factor(s) during an MD run
!> \param girsanov_env ...
!> \param itimes ...
!> \param time ...
!> \param output_unit ...
!> \param new_file ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE write_reweighting_factor(girsanov_env, itimes, time, output_unit, new_file)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      INTEGER, INTENT(IN)                                :: itimes
      REAL(KIND=dp), INTENT(IN)                          :: time
      INTEGER, INTENT(IN)                                :: output_unit
      LOGICAL, INTENT(IN)                                :: new_file

      IF (new_file) THEN
         WRITE (UNIT=output_unit, FMT='(A,22X,A,22X,A))') "#   Step   Time [fs]", "U/RT", "ΔS/D"
      END IF
      WRITE (UNIT=output_unit, FMT='(I8,F12.3,2(1X,ES25.16E2))') &
         itimes, time, girsanov_env%U/girsanov_env%kBT, girsanov_env%delta_S

      CALL m_flush(output_unit)

      ! reset S after output
      girsanov_env%delta_S = 0.0_dp

   END SUBROUTINE write_reweighting_factor

END MODULE girsanov_utils
