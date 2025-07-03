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
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type  
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_constants,                 ONLY: grw_bias,&
                                              grw_perturbation,&
                                              grw_bias_metadynamics,&
                                              grw_bias_external,&
                                              grw_bias_restraint,&
                                              grw_bias_mixing,&
                                              grw_perturbation_linres                                              
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE memory_utilities,                ONLY: reallocate
   USE particle_list_types,             ONLY: particle_list_type
   USE physcon,                         ONLY: boltzmann, &
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
                                                            temperature, dt
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

      ! CALL cite_reference(VandeVondele2002)
      NULLIFY (girsanov_section, mode_section, i_atoms, atom_list)

      check = .NOT. ASSOCIATED(girsanov_env)
      CPASSERT(check)

      girsanov_section => section_vals_get_subs_vals(md_section, "GIRSANOV")
      CALL section_vals_get(girsanov_section, explicit=explicit)

      IF (explicit) THEN
         CALL timeset(routineN, handle)
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


!      INSERT ALL THE CHECKS HERE, such as: (since we have the md_section, why not?)
!          check = barostat%virial_components == do_clv_xyz .OR. simpar%ensemble == npt_f_ensemble
!          IF (.NOT. check) &
!             CALL cp_abort(__LOCATION__, "The screening of the components of "// &
!                           "the virial is available only with the NPT_F ensemble!")
         ! - ENSEMBLE (NPT or NVT)
         ! - THEROMSTAT (CSVR MASSIVE) or at least one thermo per biased DOF (atomlist)
         ! - if MTD: PLUMED, no extended Lagrangian (if we exlude NON-PLUMED, do we exclude this was well?)

         ! DEALLOCATE(girsanov_section)
         ! DEALLOCATE(mode_section)
         CALL timestop(handle)
      END IF


   END SUBROUTINE create_girsanov_environment

! **************************************************************************************************
!> \brief initializes the Girsanov reweighting environment
!> \param girsanov_env ...
!> \par History
!>      07.2025 created [sjaehnigen]
!> \author Sascha Jähnigen - Freie Universität Berlin
! **************************************************************************************************   
   SUBROUTINE initialize_girsanov_environment(girsanov_env)
      TYPE(girsanov_env_type), POINTER             :: girsanov_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'initialize_girsanov_environment'
      INTEGER                                            :: handle    
      
      REAL(KIND=dp)                                      :: f_prime, d_prime

      IF (ASSOCIATED(girsanov_env)) THEN
         CALL timeset(routineN, handle)      

         ! calculate parameters for update operators according to OBABO scheme
         d_prime = EXP(-0.5_dp * girsanov_env%xi * girsanov_env%dt)
         f_prime = SQRT(boltzmann/joule * girsanov_env%temperature*kelvin * (1.0_dp - d_prime**2))
         ! NB: division by SQRT(mass) on the fly in girsanov_utils.F
         girsanov_env%prefactor_eta_1 = -0.5_dp * girsanov_env%dt / f_prime
         girsanov_env%prefactor_eta_2 = girsanov_env%prefactor_eta_1 * d_prime
   
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
!> \param i_step ...
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

      INTEGER                                            :: idx
      INTEGER, POINTER, DIMENSION(:, :)                  :: index_map

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(girsanov_env%degrees_of_freedom))

      dof => girsanov_env%degrees_of_freedom
      index_map => girsanov_env%index_map

      idx = index_map(i_dim, i_particle)
      IF (idx /= 0) THEN
         dof(idx)%grad_U(2) = dof(idx)%grad_U(2) + grad_U
      END IF
 
      NULLIFY(dof)
      NULLIFY(index_map)

      CALL timestop(handle)

      ! DO WE NEED MP BROADCAST HERE?
      ! my_overwrite = .FALSE.
      ! IF (PRESENT(overwrite)) my_overwrite = overwrite
   END SUBROUTINE update_girsanov_environment_r0

! **************************************************************************************************
!> \brief updates the Girsanov environment (rank-1 input for particle)
!> \param girsanov_env ...
!> \param i_particle ...
!> \param i_step ...
!> \param grad_U ... (add to what is there)
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

      INTEGER                                            :: i, j, n_dim, n_atoms, idx1, idx2
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
            idx1 = index_map(j, index)
            IF (idx1 /= 0) THEN
               dof(idx1)%grad_U(2) = dof(idx1)%grad_U(2) + grad_U(j)
            END IF
         END DO
      ELSE
         DO i=1, n_atoms
            idx1 = (i-1)*n_dim+index
            idx2 = atom_list(i)
            dof(idx1)%grad_U(2) = dof(idx1)%grad_U(2) + grad_U(idx2)
         END DO         
      END IF

      NULLIFY(dof)
      NULLIFY(index_map)
      NULLIFY(atom_list)

      CALL timestop(handle)

      ! DO WE NEED MP BROADCAST HERE?
      ! my_overwrite = .FALSE.
      ! IF (PRESENT(overwrite)) my_overwrite = overwrite
   END SUBROUTINE update_girsanov_environment_r1

! **************************************************************************************************
!> \brief updates the Girsanov environment (rank-2 input)
!> \param girsanov_env ...
!> \param i_step ...
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

      INTEGER                                            :: i, j, n_atoms, n_dim, idx
      INTEGER, POINTER, DIMENSION(:)                     :: atom_list

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(girsanov_env%degrees_of_freedom))

      dof => girsanov_env%degrees_of_freedom
      n_atoms = girsanov_env%n_atoms
      n_dim = girsanov_env%n_dim
      atom_list => girsanov_env%atom_list

      DO j=1, n_dim
         DO i=1, n_atoms
            idx = atom_list(i)
            dof((i-1)*n_dim+j)%grad_U(2) = dof((i-1)*n_dim+j)%grad_U(2) + grad_U(j, idx)
         END DO
      END DO 

      NULLIFY(dof)
      NULLIFY(atom_list)

      CALL timestop(handle)

      ! DO WE NEED MP BROADCAST HERE?
      ! my_overwrite = .FALSE.
      ! IF (PRESENT(overwrite)) my_overwrite = overwrite
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

END MODULE girsanov_types
