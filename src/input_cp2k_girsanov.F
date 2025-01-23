!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input section for Girsanov reweighting
!> \par History
!>      09.2024 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
MODULE input_cp2k_girsanov
!      USE bibliography,                    ONLY: &
!         Aguado2003, BaniHashemian2016, Blochl1995, Darden1993, Essmann1995, Ewald1921, &
!         Genovese2006, Genovese2007, Laino2008, Martyna1999, Toukmaji1996

!    USE barostat_types,                  ONLY: do_clv_x,&
!                                               do_clv_xy,&
!                                               do_clv_xyz,&
!                                               do_clv_xz,&
!                                               do_clv_y,&
!                                               do_clv_yz,&
!                                               do_clv_z
   USE input_constants,                 ONLY: grw_bias,&
                                              grw_perturbation,&
                                              grw_bias_use_metadynamics,&
                                              grw_bias_use_external,&
                                              grw_bias_use_mixing,&
                                              grw_perturbation_use_linres
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              high_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
!    USE input_cp2k_thermostats,          ONLY: create_mass_section,&
!                                               create_thermostat_section,&
!                                               create_velocity_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
!    USE input_val_types,                 ONLY: real_t
!    USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: newline, s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_girsanov'

   PUBLIC :: create_girsanov_section

!***
CONTAINS

! **************************************************************************************************
!> \brief Create Girsanov section
!> \param section ...
!> \author Sascha Jähnigen - Freie Universität Berlin - 09.2024
! **************************************************************************************************
   SUBROUTINE create_girsanov_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: bias_section, perturbation_section, subsection
      ! subsection, thermo_section

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="GIRSANOV", &
                          description="Parameters for Girsanov reweighting."//newline//"Hier könnte Ihre Werbung stehen.", &
                          ! citations=(/Kieninger2023, Schäfer2024/), &
                          n_keywords=1, n_subsections=3, repeats=.FALSE.)

      NULLIFY (keyword, bias_section, perturbation_section, subsection)
      CALL keyword_create(keyword, __LOCATION__, name="MODE", &
                          description="Source of the bias potential U", &
                          usage="MODE (BIAS | PERTURBATION)", &
                          enum_c_vals=s2a("BIAS", "PERTURBATION"), &
                          enum_i_vals=(/grw_bias, grw_perturbation/), &
                          enum_desc=s2a("bias", "perturbation"), &
                          default_i_val=grw_bias)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_bias_section(bias_section)
      CALL section_add_subsection(section, bias_section)
      CALL section_release(bias_section)

      CALL create_perturbation_section(perturbation_section)
      CALL section_add_subsection(section, perturbation_section)
      CALL section_release(perturbation_section)

      CALL create_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_girsanov_section

! **************************************************************************************************
!> \brief Create bias subsection for Girsanov section
!> \param section ...
!> \author Sascha Jähnigen - Freie Universität Berlin - 12.2024
! **************************************************************************************************
   SUBROUTINE create_bias_section(section)
      TYPE(section_type), POINTER                        :: section
      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="BIAS", &
                          description="Parameters of Girsanov reweighting using a bias potential that has been applied to "// &
                          "the system. This section requires that the corresponding source method be enabled.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOURCE", &
                          description="Specify the source from which to extract the bias potential", &
                          usage="SOURCE (METADYNAMICS | EXTERNAL_POTENTIAL | MIXING)", &
                          enum_c_vals=s2a("METADYNAMICS", "EXTERNAL_POTENTIAL", "MIXING"), &
                          enum_i_vals=(/grw_bias_use_metadynamics, grw_bias_use_external, grw_bias_use_mixing/), &
                          enum_desc=s2a("Use hills from metadynamics. Requires the METADYN section etc.", &
                                        "Use external potential. Requires the EXTERNAL_POTENTIAL section.", &
                                        "Use potential from force mixing. Requires the MIXED section etc."), &
                          default_i_val=grw_bias_use_metadynamics)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_bias_section

! **************************************************************************************************
!> \brief Create perturbation subsection for Girsanov section
!> \param section ...
!> \author Sascha Jähnigen - Freie Universität Berlin - 12.2024
! **************************************************************************************************
   SUBROUTINE create_perturbation_section(section)
      TYPE(section_type), POINTER                        :: section
      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PERTURBATION", &
                          description="Parameters of Girsanov reweighting using a perturbation potential that has NOT been "// &
                          "applied to the system. This section requires that the corresponding source method be enabled.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SOURCE", &
                          description="Specify the source of the perturbation potential", &
                          usage="SOURCE (LINRES)", &
                          enum_c_vals=s2a("LINRES"), &
                          enum_i_vals=(/grw_perturbation_use_linres/), &
                          enum_desc=s2a("Use results from a linear response calculation. Requires the LINRES section."), &
                          default_i_val=grw_perturbation_use_linres)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_perturbation_section

! **************************************************************************************************
!> \brief Create print subsection for Girsanov section
!> \param section ...
!> \author Sascha Jähnigen - Freie Universität Berlin - 09.2024
! **************************************************************************************************
   SUBROUTINE create_print_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (print_key)
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Collects all print keys for Girsanov reweighting.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "REWEIGHTING_FACTORS", &
                                       description="Controls the output of the static and dynamic reweighting factors "// &
                                       "for all biased degrees of freedom.", &
                                       print_level=high_print_level, common_iter_levels=1, &
                                       filename="")
      ! Switch bewteen RNG and RW factors through format keyword? (example taken from velocities)
      ! CALL add_format_keyword(keyword, print_key, pos=.FALSE., &
      !                         description="Specifies the format of the output file for the velocities.")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)
   END SUBROUTINE create_print_section

END MODULE input_cp2k_girsanov
