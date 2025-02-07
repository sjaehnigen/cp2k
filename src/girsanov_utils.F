!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
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
   USE girsanov_types,                  ONLY: girsanov_env_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_constants,                 ONLY: grw_bias,&
                                              grw_perturbation
   USE kinds,                           ONLY: dp
                                              !default_string_length,&
                                              !dp
!   USE machine,                         ONLY: m_flush
!   USE physcon,                         ONLY: angstrom,&
!                                              femtoseconds,&
!                                              kelvin
!   USE simpar_types,                    ONLY: simpar_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC ::  write_reweighting_factor

! *** Global parameters ***
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'girsanov_utils'

CONTAINS
! **************************************************************************************************
!> \brief Write Girsanov reweighting factor(s) during an MD run
!> \param barostat ...
!> \param simpar ...
!> \param my_pos ...
!> \param my_act ...
!> \param cell ...
!> \param itimes ...
!> \param time ...
!> \par History
!>      12.2024 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE write_reweighting_factor(girsanov_env) ! barostat, simpar, my_pos, my_act, cell, itimes, time)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
!      TYPE(barostat_type), POINTER                       :: barostat
!      TYPE(simpar_type), INTENT(IN)                      :: simpar
!      CHARACTER(LEN=default_string_length)               :: my_pos, my_act
!      TYPE(cell_type), POINTER                           :: cell
!      INTEGER, INTENT(IN)                                :: itimes
      ! REAL(KIND=dp), INTENT(IN)                          :: time
!
!      INTEGER                                            :: baro, nfree
!      LOGICAL                                            :: new_file
!      REAL(KIND=dp)                                      :: baro_kin, baro_pot, temp
!      TYPE(cp_logger_type), POINTER                      :: logger
!
!      NULLIFY (logger)
!      logger => cp_get_default_logger()
!
!      IF (ASSOCIATED(barostat)) THEN
!         baro = cp_print_key_unit_nr(logger, barostat%section, "PRINT%ENERGY", &
!                                     extension=".bener", file_position=my_pos, file_action=my_act, is_new_file=new_file)
!         CALL get_baro_energies(cell, simpar, barostat%npt, baro_kin, baro_pot)
!         nfree = SIZE(barostat%npt, 1)*SIZE(barostat%npt, 2)
!         temp = 2.0_dp*baro_kin/REAL(nfree, dp)*kelvin
!         IF (baro > 0) THEN
!            IF (new_file) THEN
!               WRITE (baro, '("#",3X,A,10X,A,8X,3(5X,A,5X),3X,A)') "Step Nr.", "Time[fs]", "Kin.[a.u.]", &
!                  "Temp.[K]", "Pot.[a.u.]", "Vol[Ang.^3]"
!            END IF
!            WRITE (UNIT=baro, FMT="(I10, F20.3,4F20.10)") itimes, time*femtoseconds, &
!               baro_kin, temp, baro_pot, cell%deth*angstrom*angstrom*angstrom
!            CALL m_flush(baro)
!         END IF
!         CALL cp_print_key_finished_output(baro, logger, barostat%section, "PRINT%ENERGY")
!      END IF
      IF (ASSOCIATED(girsanov_env)) THEN
          CALL cp_warn(__LOCATION__, &
                     ">>> I was told to write Girsanov output! "// &
                     "Congratulations!")
      END IF

   END SUBROUTINE write_reweighting_factor

END MODULE girsanov_utils
