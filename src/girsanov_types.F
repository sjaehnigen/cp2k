!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type  
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type 
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr                                                
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_constants,                 ONLY: grw_bias,&
                                              grw_perturbation,&
                                              grw_bias_metadynamics,&
                                              grw_bias_external,&
                                              grw_bias_restraint,&
                                              grw_bias_mixing,&
                                              grw_perturbation_linres,&
                                              nvt_ensemble                                            
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: boltzmann, &
                                              femtoseconds, &
                                              joule, &
                                              kelvin
   USE simpar_types,                    ONLY: simpar_type

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: girsanov_env_type, &
             degree_of_freedom_type, &
             create_girsanov_environment, &
             initialize_girsanov_environment, &
             roll_girsanov_environment, &
             release_girsanov_environment, &
             update_girsanov_environment

   INTERFACE update_girsanov_environment
      MODULE PROCEDURE update_girsanov_environment_static, &
                       update_girsanov_environment_r0, &
                       update_girsanov_environment_r1, &
                       update_girsanov_environment_r2
   END INTERFACE update_girsanov_environment       

! **************************************************************************************************
!> \brief defines the data type for calculating for a single DOF the reweighting factor and buffer 
! **************************************************************************************************
   TYPE degree_of_freedom_type
      ! grad U of the last two time steps (2->last, 1->previous), i.e. before and after A step
      REAL(KIND=dp), DIMENSION(2)                        :: grad_U = 0.0_dp
      ! random numbers of the two O steps
      REAL(KIND=dp), DIMENSION(2)                        :: eta = 0.0_dp 
      REAL(KIND=dp), DIMENSION(2)                        :: delta_eta = 0.0_dp 
      REAL(KIND=dp)                                      :: mass
   END TYPE 

! **************************************************************************************************
!> \brief defines the environment type for Girsanov reweighting
!> \param atom_list ... list of atoms for which three reweighting factors are calculated (1 per DOF)
!> \param index_map ... maps dimension and particle_index to internal running index of girsanov_env
!> \param thermo_map ... maps thermostat_index (=thermostatted DOF) to running index of girsanov_env
! **************************************************************************************************
   TYPE girsanov_env_type
      INTEGER                                            :: mode, n_atoms, n_dim
      INTEGER, POINTER, DIMENSION(:)                     :: atom_list => NULL()
      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map => NULL()
      INTEGER, POINTER, DIMENSION(:)                     :: thermo_map => NULL()
      REAL(KIND=dp)                                      :: U = 0.0_dp
      REAL(KIND=dp)                                      :: delta_S = 0.0_dp
      REAL(KIND=dp)                                      :: xi, prefactor_eta_1, prefactor_eta_2,&
                                                            temperature, dt, kBT
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: degrees_of_freedom => NULL()
   END TYPE

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
   SUBROUTINE create_girsanov_environment(girsanov_env, md_section, simpar, subsys, n_dim)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(section_vals_type), POINTER                   :: md_section
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, OPTIONAL                                  :: n_dim
      INTEGER                                            :: grw_mode, grw_source, &
                                                            n_rep, n_atoms, n_particles, &
                                                            k, i, j, n
      INTEGER, DIMENSION(:), POINTER                     :: i_atoms, atom_list
      
      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_girsanov_environment'
      INTEGER                                            :: handle      

      LOGICAL                                            :: check, explicit, index_in_range
      TYPE(section_vals_type), POINTER                   :: girsanov_section, mode_section

      ! CALL cite_reference(our papers)
      NULLIFY (girsanov_section, mode_section, i_atoms, atom_list)

      check = .NOT. ASSOCIATED(girsanov_env)
      CPASSERT(check)

      girsanov_section => section_vals_get_subs_vals(md_section, "GIRSANOV")
      CALL section_vals_get(girsanov_section, explicit=explicit)

      IF (explicit) THEN
         CALL timeset(routineN, handle)

         IF (simpar%ensemble /= nvt_ensemble) THEN
            CPABORT("Girsanov reweighting not implemented for this ensemble (use NVT).")
         END IF

         ALLOCATE (girsanov_env)
         girsanov_env%dt = simpar%dt
         girsanov_env%temperature = simpar%temp_ext

         IF (PRESENT(n_dim)) THEN
            girsanov_env%n_dim = n_dim
         ELSE
            girsanov_env%n_dim = 3
         END IF                     
         CALL section_vals_val_get(girsanov_section, "MODE", i_val=grw_mode)
         SELECT CASE (grw_mode)
         CASE (grw_bias)
            mode_section => section_vals_get_subs_vals(girsanov_section, "BIAS")
         CASE (grw_perturbation)
            mode_section => section_vals_get_subs_vals(girsanov_section, "PERTURBATION")   
         CASE DEFAULT
            CPABORT("Unknown MODE in GIRSANOV section")
         END SELECT           

         CALL section_vals_val_get(mode_section, "SOURCE", i_val=grw_source)
         SELECT CASE (grw_source)
         CASE(grw_bias_metadynamics)
            girsanov_env%mode = grw_bias_metadynamics
         CASE(grw_bias_external)
            girsanov_env%mode = grw_bias_external
         CASE(grw_bias_restraint)
            girsanov_env%mode = grw_bias_restraint
         CASE(grw_bias_mixing)
            CPABORT("Girsanov reweighting not yet implemented for bias from mixed force environments")
         CASE(grw_perturbation_linres)
            CPABORT("Girsanov reweighting not yet implemented for perturbation from linear response")
         CASE DEFAULT
            CPABORT("Unknown SOURCE in GIRSANOV section")
         END SELECT    
         
         ! ATOM LIST 
         n_atoms = 0
         CALL section_vals_val_get(mode_section, "ATOM_LIST", n_rep_val=n_rep)
         DO k = 1, n_rep
            CALL section_vals_val_get(mode_section, "ATOM_LIST", i_rep_val=k, i_vals=i_atoms)
            CALL reallocate(atom_list, 1, n_atoms + SIZE(i_atoms))
            atom_list(n_atoms + 1:n_atoms + SIZE(i_atoms)) = i_atoms
            n_atoms = n_atoms + SIZE(i_atoms)
         END DO

         CALL cp_subsys_get(subsys, &
                            particles=particles)
         n_particles = particles%n_els
         ALLOCATE(girsanov_env%index_map(girsanov_env%n_dim, n_particles))
         girsanov_env%index_map = 0

         ! use all atoms
         IF ((n_rep == 0) .OR. ANY(atom_list == -1)) THEN
            n_atoms = n_particles
            CALL reallocate(atom_list, 1, n_atoms)
            atom_list = [(i, i = 1, n_atoms)]
            DO i = 1, n_atoms
               ! create (dummy) hash table for inverse atom list search
               DO n = 1, girsanov_env%n_dim
                  girsanov_env%index_map(n, i) = (i-1)*girsanov_env%n_dim + n
               END DO
            END DO            
         ELSE
            ! clean and verify atom_list
            CPASSERT(SIZE(atom_list) /= 0)
            index_in_range = (MAXVAL(atom_list) <= n_particles) .AND. (MINVAL(atom_list) > 0)
            IF (.NOT. index_in_range) &
               CPABORT("Index out of range in atom list for GIRSANOV reweighting.")
            DO i = 1, n_atoms
               DO j = i + 1, n_atoms
                  IF (atom_list(i) == atom_list(j)) &
                     CPABORT("Found duplicate entries in atom list for GIRSANOV reweighting.")
               END DO
               ! create hash table for inverse atom list search
               DO n=1, girsanov_env%n_dim
                  girsanov_env%index_map(n, atom_list(i)) = (i-1)*girsanov_env%n_dim + n
               END DO
            END DO
         END IF
         
         ! set up lists
         ALLOCATE(girsanov_env%degrees_of_freedom(girsanov_env%n_dim*n_atoms))
         girsanov_env%atom_list => atom_list
         girsanov_env%n_atoms = n_atoms
         
         ! fill with particle info
         DO i = 1, girsanov_env%n_atoms
            DO  n = 1, girsanov_env%n_dim
               j = girsanov_env%index_map(n, atom_list(i))
               girsanov_env%degrees_of_freedom(j)%mass = particles%els(atom_list(i))%atomic_kind%mass
            END DO
         ENDDO

         CALL timestop(handle)
      END IF


   END SUBROUTINE create_girsanov_environment

! **************************************************************************************************
!> \brief initializes the Girsanov reweighting environment
!> \param girsanov_env ...
!> \param md_section ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************   
   SUBROUTINE initialize_girsanov_environment(girsanov_env, md_section)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(section_vals_type), POINTER                   :: md_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'initialize_girsanov_environment'
      INTEGER                                            :: handle    
      
      TYPE(section_vals_type), POINTER                   :: girsanov_section
      REAL(KIND=dp)                                      :: f_prime, d_prime

      IF (ASSOCIATED(girsanov_env)) THEN
         CALL timeset(routineN, handle)      

         girsanov_section => section_vals_get_subs_vals(md_section, "GIRSANOV")

         girsanov_env%kBT = boltzmann/joule * girsanov_env%temperature*kelvin
         ! calculate parameters for update operators according to OBABO scheme
         d_prime = EXP(-0.5_dp * girsanov_env%xi * girsanov_env%dt)
         f_prime = SQRT(girsanov_env%kBT * (1.0_dp - d_prime**2))
         ! NB: division by SQRT(mass) on the fly in girsanov_utils.F
         girsanov_env%prefactor_eta_1 = 0.5_dp * girsanov_env%dt / f_prime
         girsanov_env%prefactor_eta_2 = girsanov_env%prefactor_eta_1 * d_prime
   
         CALL girsanov_info(girsanov_env, girsanov_section)

         ! add print info to outfile
         CALL timestop(handle)
      END IF
      
   END SUBROUTINE initialize_girsanov_environment

! **************************************************************************************************
!> \brief refreshes DOF array: moves grad_U: 2->1 and resets eta and U
!> \param girsanov_env ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************   
   SUBROUTINE roll_girsanov_environment(girsanov_env)
      TYPE(girsanov_env_type), INTENT(INOUT), POINTER    :: girsanov_env
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: dof
      ! INTEGER                                            :: i

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'roll_girsanov_environment'
      INTEGER                                            :: handle      

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(girsanov_env))

      dof => girsanov_env%degrees_of_freedom

      dof(:)%grad_U(1) = dof(:)%grad_U(2)
      dof(:)%grad_U(2) = 0.0_dp
      dof(:)%eta(1) = 0.0_dp
      dof(:)%eta(2) = 0.0_dp
      girsanov_env%U = 0.0_dp
      NULLIFY(dof)
            
      CALL timestop(handle)

   END SUBROUTINE roll_girsanov_environment

! **************************************************************************************************
!> \brief updates the Girsanov environment (static input)
!> \param girsanov_env ...
!> \param U ... (add to what is there)
!> \par History
!>      06.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE update_girsanov_environment_static(girsanov_env, U)   
      TYPE(girsanov_env_type), INTENT(INOUT), POINTER    :: girsanov_env
      REAL(KIND=dp), INTENT(IN)                          :: U

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_girsanov_environment_static'
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(girsanov_env))

      girsanov_env%U = girsanov_env%U + U
      
      CALL timestop(handle)

   END SUBROUTINE update_girsanov_environment_static


! **************************************************************************************************
!> \brief updates the Girsanov environment (rank-0 input)
!> \param girsanov_env ...
!> \param i_particle ...
!> \param i_dim ...
!> \param grad_U ... (add to what is there)
!> \par History
!>      02.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE update_girsanov_environment_r0(girsanov_env, i_particle, i_dim, grad_U)   
      TYPE(girsanov_env_type), INTENT(INOUT)             :: girsanov_env
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: dof
      INTEGER, INTENT(IN)                                :: i_particle, i_dim
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: grad_U

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_girsanov_environment_r0'
      INTEGER                                            :: handle

      INTEGER                                            :: i_entry
      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(girsanov_env%degrees_of_freedom))

      dof => girsanov_env%degrees_of_freedom
      index_map => girsanov_env%index_map

      i_entry = index_map(i_dim, i_particle)
      IF (i_entry /= 0) THEN
         dof(i_entry)%grad_U(2) = dof(i_entry)%grad_U(2) + grad_U
      END IF
 
      NULLIFY(dof)
      NULLIFY(index_map)

      CALL timestop(handle)

   END SUBROUTINE update_girsanov_environment_r0

! **************************************************************************************************
!> \brief updates the Girsanov environment (rank-1 input for particle or dimension)
!> \param girsanov_env ...
!> \param index ...
!> \param grad_U ... (add to what is there)
!> \param index_is_particle
!> \par History
!>      02.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE update_girsanov_environment_r1(girsanov_env, index, grad_U, index_is_particle) 
      TYPE(girsanov_env_type), INTENT(INOUT)             :: girsanov_env
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER  &
                                                         :: dof
      INTEGER, INTENT(IN)                                :: index
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: grad_U
      LOGICAL, INTENT(IN), OPTIONAL                      :: index_is_particle

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_girsanov_environment_r1'
      INTEGER                                            :: handle

      INTEGER                                            :: i, j, n_dim, n_atoms, i_entry, i_particle
      INTEGER, POINTER, DIMENSION(:)                     :: atom_list
      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map
      LOGICAL                                            :: index_is_dim

      CALL timeset(routineN, handle)
      
      CPASSERT(ASSOCIATED(girsanov_env%degrees_of_freedom))

      dof => girsanov_env%degrees_of_freedom
      n_dim = girsanov_env%n_dim
      n_atoms = girsanov_env%n_atoms
      index_map => girsanov_env%index_map
      atom_list => girsanov_env%atom_list
      
      index_is_dim = .True.
      IF (PRESENT(index_is_particle) .AND. index_is_particle) &
         index_is_dim = .False.

      IF (.NOT. index_is_dim) THEN
         DO j=1, n_dim
            i_entry = index_map(j, index)
            IF (i_entry /= 0) THEN
               dof(i_entry)%grad_U(2) = dof(i_entry)%grad_U(2) + grad_U(j)
            END IF
         END DO
      ELSE
         DO i=1, n_atoms
            i_particle = atom_list(i)
            i_entry = index_map(index, i_particle)   
            dof(i_entry)%grad_U(2) = dof(i_entry)%grad_U(2) + grad_U(i_particle)
         END DO         
      END IF

      NULLIFY(dof)
      NULLIFY(index_map)
      NULLIFY(atom_list)

      CALL timestop(handle)

   END SUBROUTINE update_girsanov_environment_r1

! **************************************************************************************************
!> \brief updates the Girsanov environment (rank-2 input)
!> \param girsanov_env ...
!> \param grad_U ... (add to what is there)
!> \par History
!>      02.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE update_girsanov_environment_r2(girsanov_env, grad_U)   
      TYPE(girsanov_env_type), INTENT(INOUT)             :: girsanov_env
      TYPE(degree_of_freedom_type), DIMENSION(:), POINTER &
                                                         :: dof
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), OPTIONAL &
                                                         :: grad_U

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'update_girsanov_environment_r2'
      INTEGER                                            :: handle

      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map
      INTEGER                                            :: i, j, n_atoms, n_dim, i_entry, i_particle
      INTEGER, POINTER, DIMENSION(:)                     :: atom_list

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(girsanov_env%degrees_of_freedom))

      dof => girsanov_env%degrees_of_freedom
      n_atoms = girsanov_env%n_atoms
      index_map => girsanov_env%index_map
      n_dim = girsanov_env%n_dim
      atom_list => girsanov_env%atom_list

      DO j=1, n_dim
         DO i=1, n_atoms
            i_particle = atom_list(i)
            i_entry = index_map(j, i_particle) 
            dof(i_entry)%grad_U(2) = dof(i_entry)%grad_U(2) + grad_U(j, i_particle)
         END DO
      END DO 

      NULLIFY(dof)
      NULLIFY(index_map)
      NULLIFY(atom_list)

      CALL timestop(handle)

   END SUBROUTINE update_girsanov_environment_r2
   
! **************************************************************************************************
!> \brief releases the Girsanov environment
!> \param girsanov_env ...
!> \par History
!>      01.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE release_girsanov_environment(girsanov_env)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env

      IF (ASSOCIATED(girsanov_env)) THEN
          ! NULLIFY SECTION
         IF (ASSOCIATED(girsanov_env%atom_list)) DEALLOCATE(girsanov_env%atom_list)
         IF (ASSOCIATED(girsanov_env%index_map)) DEALLOCATE(girsanov_env%index_map)
         IF (ASSOCIATED(girsanov_env%thermo_map)) DEALLOCATE(girsanov_env%thermo_map)
         IF (ASSOCIATED(girsanov_env%degrees_of_freedom)) DEALLOCATE(girsanov_env%degrees_of_freedom)
         DEALLOCATE (girsanov_env)
      END IF

   END SUBROUTINE release_girsanov_environment

! **************************************************************************************************
!> \brief prints info on Girsanov environment
!> \param girsanov_env ...
!> \param girsanov_section ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************
   SUBROUTINE girsanov_info(girsanov_env, girsanov_section)
      TYPE(girsanov_env_type), POINTER                   :: girsanov_env
      TYPE(section_vals_type), POINTER                   :: girsanov_section

      INTEGER                                            :: iw
      TYPE(cp_logger_type), POINTER                      :: logger
      
      TYPE(girsanov_env_type), POINTER                   :: grw
      CHARACTER(len=18)                                  :: token, mode, source
      CHARACTER(len=1000)                                :: output
      INTEGER                                            :: i, start_idx, end_idx, linelen
      INTEGER, PARAMETER                                 :: max_line_len = 30

      NULLIFY (logger, grw)
      logger => cp_get_default_logger()
      iw = cp_print_key_unit_nr(logger, girsanov_section,&
                                "PRINT%GIRSANOV_INFO", extension=".girsanov_log")
      
      grw => girsanov_env
      IF (iw > 0) THEN
         WRITE (iw, '(/,T2,A)') 'GIRSANOV| Dynamical reweighting information'
         SELECT CASE (grw%mode)
         CASE(grw_bias_metadynamics)
            mode = "Bias"
            source = "Metadynamics"
         CASE(grw_bias_external)
            mode = "Bias"
            source = "External Potential"
         CASE(grw_bias_restraint)
            mode = "Bias"
            source = "Restraint"
         END SELECT
         WRITE (iw, '(T2,A,T57,A24)') 'GIRSANOV| Source ', TRIM(mode)//": "//TRIM(source)
         WRITE (iw, '(T2,A,T63,F18.6)') 'GIRSANOV| Gamma [fs^-1]', grw%xi/femtoseconds
         ! WRITE (iw, '(T2,A,T63,F18.1)') 'GIRSANOV| Temperature ', grw%temperature*kelvin
         WRITE (iw, '(T2,A,T71,I10)') 'GIRSANOV| Number of monitored degrees of freedom ', SIZE(grw%degrees_of_freedom)
         i = 1
         output = ""
         linelen = 0
         DO WHILE (i <= grw%n_atoms)
            start_idx = i
            DO WHILE (i < grw%n_atoms .AND. grw%atom_list(i+1) == grw%atom_list(i) + 1)
               i = i + 1
            END DO
            end_idx = i

            IF (end_idx > start_idx + 1) THEN
               WRITE(token, '(I0,"..",I0)') grw%atom_list(start_idx), grw%atom_list(end_idx)
            ELSE IF (end_idx == start_idx + 1) THEN
               WRITE(token, '(I0," ",I0)') grw%atom_list(start_idx), grw%atom_list(end_idx)
            ELSE
               WRITE(token, '(I0)') grw%atom_list(start_idx)
            END IF

            IF (linelen > 0) token = ' ' // TRIM(token)
            IF (linelen + LEN_TRIM(token) > max_line_len) THEN
               WRITE (iw, '(T2,A,T51,A30)') 'GIRSANOV| Atom list ', TRIM(output)
               output = TRIM(token)
            ELSE
               output = TRIM(output) // TRIM(token)
            END IF
            linelen = LEN_TRIM(output)

            i = i + 1
         END DO
         WRITE (iw, '(T2,A,T51,A30)') 'GIRSANOV| Atom list ', TRIM(output)
         
      END IF

      CALL cp_print_key_finished_output(iw, logger, girsanov_section, "PRINT%GIRSANOV_INFO")

   END SUBROUTINE girsanov_info

END MODULE girsanov_types
