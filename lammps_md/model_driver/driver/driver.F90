!
! CDDL HEADER START
!
! The contents of this file are subject to the terms of the Common Development
! and Distribution License Version 1.0 (the "License").
!
! You can obtain a copy of the license at
! http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
! specific language governing permissions and limitations under the License.
!
! When distributing Covered Code, include this CDDL HEADER in each file and
! include the License file in a prominent location with the name LICENSE.CDDL.
! If applicable, add the following below this CDDL HEADER, with the fields
! enclosed by brackets "[]" replaced with your own identifying information:
!
! Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
!
! CDDL HEADER END
!
!
! Copyright (c) 2013--2018.
! All rights reserved.
!
! Contributors:
!    Alireza Khorshidi
!    C. Franklin Goldsmith
!    Ryan S. Elliott
!    Malte Doentgen
!    Muammar El-Khatib Rodriguez
!
module my_driver

use, intrinsic :: iso_c_binding

use char_zeta
use neuralnetwork
use kim_model_driver_headers_module

implicit none

save
private
public BUFFER_TYPE,               &
       Compute_Energy_Forces,     &
       compute_arguments_create,  &
       compute_arguments_destroy, &
       refresh,                   &
       destroy
integer(c_int), parameter :: cd = c_double 
integer(c_int), parameter :: DIM=3

type:: integer_one_d_array
  sequence
  integer(c_int), allocatable:: onedarray(:)
end type integer_one_d_array

type :: atomic_descriptors
  character(len=20, kind=c_char) :: symmetry_type
  integer(c_int) :: atomic_num_first, atomic_num_second, atomic_num_third
  real(c_double) :: g_gamma, g_zeta
  integer :: num_rscalars, num_rtypes
  real(c_double), allocatable :: r_scalars(:)
  integer, allocatable :: r_types(:)
end type atomic_descriptors

type :: neural_network
  integer(c_int) :: no_layers
  integer(c_int), allocatable :: no_nodes(:)
  real(c_double), allocatable :: parameters(:)
end type neural_network

type :: model
  ! navigator --> which NN should be used for an element type
  integer(c_int), allocatable :: navigator(:)
  type(neural_network), allocatable ::  networks(:)
  real(c_double), allocatable :: min_fingerprints(:)
  real(c_double), allocatable :: max_fingerprints(:)
end type model

!-------------------------------------------------------------------------------
!  Definition of Buffer type
!-------------------------------------------------------------------------------
type :: BUFFER_TYPE
  real(c_double) :: influence_distance(1)
  real(c_double) :: cutoff(1)
  integer(c_int) :: not_requesting_n_from_ghost(1)
  
  integer(c_int) :: num_species
  character(len=8, kind=c_char), allocatable :: species_name_strings(:)
  type(kim_species_name_type), allocatable :: species_names(:)
  integer(c_int), allocatable :: atomic_numbers(:)
  integer(c_int), allocatable :: species_codes(:)
  
  integer(c_int):: cutofffn_code, p_gamma
  integer(c_int) :: num_descriptors
  ! one block corresponds to 
  ! for a fixed tuple of (gamma, zeta), a list of r_scalars and r_types
  integer(c_int) :: num_blocks
  type(atomic_descriptors), allocatable :: descriptors(:)
  
  integer(c_int) :: ensemble_size
  type(model), allocatable :: model_ensemble(:)
end type BUFFER_TYPE

contains

!-------------------------------------------------------------------------------
! Compute energy and forces on particles from the positions.
!-------------------------------------------------------------------------------
recursive subroutine Compute_Energy_Forces(& 
  model_compute_handle, &
  model_compute_arguments_handle, ierr) bind(c)

  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_handle_type), intent(in) :: &
    model_compute_arguments_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  integer(c_int) :: ierr2
  integer(c_int) :: comp_forces,comp_energy,comp_enepot
  type(BUFFER_TYPE), pointer :: buf; type(c_ptr) :: pbuf

  !-- KIM variables
  real(c_double) :: model_cutoff
  integer(c_int), pointer :: num_atoms
  real(c_double), pointer :: energy
  real(c_double), pointer :: coor(:,:)
  real(c_double), pointer :: forces(:,:)
  real(c_double), pointer :: enepot(:)
  integer(c_int), pointer :: neighbors_of_particle(:)
  integer(c_int), pointer :: particle_species_codes(:)
  integer(c_int), pointer :: particle_contributing(:)

  integer(c_int) :: aux_ind
  integer(c_int) :: num_gs
  integer(c_int) :: num_blocks, block_count
  integer(c_int) :: cindex, csymbol
  integer(c_int) :: num_n
  integer(c_int) :: n_ind, nindex
  integer(c_int) :: force_dir

  integer(c_int) :: ensemble_count
  integer(c_int) :: model_ind

  integer(c_int), allocatable :: atomic_number_codes(:)
  integer(c_int), allocatable:: atomic_num_n(:)
  integer(c_int), allocatable:: ind_n(:)
  real(c_double), allocatable:: pos_n(:,:)
  integer(c_int) :: actual_natoms

  ! act as placeholders
  real(c_double), allocatable:: fingerprint(:), fingerprintprime(:)
  type(integer_one_d_array), allocatable :: neighborlists(:)
  real(c_double) :: atom_energy
  real(c_double):: dforce
  
  real(c_double), dimension(3):: ri
    
  ! get model buffer from KIM object
  call kim_get_model_buffer_pointer(model_compute_handle, pbuf)
  call c_f_pointer(pbuf, buf)
  
  model_cutoff = buf%influence_distance(1)
  num_gs = buf%num_descriptors
  num_blocks = buf%num_blocks
  
  ! Unpack data from KIM object
  ierr = 0
  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_number_of_particles, num_atoms, ierr2)
  ierr = ierr + ierr2
      
  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_species_codes, &
    num_atoms, particle_species_codes, ierr2)
  
  allocate(atomic_number_codes(size(particle_species_codes)))
  
! LiGePS
  do cindex = 1, size(particle_species_codes)
    if (particle_species_codes(cindex) .eq. 1) then
      atomic_number_codes(cindex) = 3
    elseif (particle_species_codes(cindex) .eq. 2) then
      atomic_number_codes(cindex) = 32
    elseif (particle_species_codes(cindex) .eq. 3) then
      atomic_number_codes(cindex) = 15
    elseif (particle_species_codes(cindex) .eq. 4) then
      atomic_number_codes(cindex) = 16
    end if
  end do

! LiPSBr
!  do cindex = 1, size(particle_species_codes)
!    if (particle_species_codes(cindex) .eq. 1) then
!      atomic_number_codes(cindex) = 3
!    elseif (particle_species_codes(cindex) .eq. 2) then
!      atomic_number_codes(cindex) = 15
!    elseif (particle_species_codes(cindex) .eq. 3) then
!      atomic_number_codes(cindex) = 16
!    elseif (particle_species_codes(cindex) .eq. 4) then
!      atomic_number_codes(cindex) = 35
!    end if
!  end do

! LiPSCl
!  do cindex = 1, size(particle_species_codes)
!    if (particle_species_codes(cindex) .eq. 1) then
!      atomic_number_codes(cindex) = 3
!    elseif (particle_species_codes(cindex) .eq. 2) then
!      atomic_number_codes(cindex) = 15
!    elseif (particle_species_codes(cindex) .eq. 3) then
!      atomic_number_codes(cindex) = 16
!    elseif (particle_species_codes(cindex) .eq. 4) then
!      atomic_number_codes(cindex) = 17
!    end if
!  end do
  
! LiN
!  do cindex = 1, size(particle_species_codes)
!    if (particle_species_codes(cindex) .eq. 1) then
!      atomic_number_codes(cindex) = 3
!    elseif (particle_species_codes(cindex) .eq. 2) then
!      atomic_number_codes(cindex) = 7
!    end if
!  end do

  ierr = ierr + ierr2

  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_particle_contributing, & 
    num_atoms, particle_contributing, &
    ierr2)
  ierr = ierr + ierr2

  actual_natoms = 0
  do cindex = 1, num_atoms
    if (particle_contributing(cindex) .eq. 1) then
      actual_natoms = actual_natoms + 1
    end if
  end do

  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_coordinates, dim, num_atoms, coor, ierr2)
  ierr = ierr + ierr2

  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_energy, energy, ierr2)
  ierr = ierr + ierr2
  
  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_forces, & 
    dim, num_atoms, forces, ierr2)
  ierr = ierr + ierr2

  call kim_get_argument_pointer( &
    model_compute_arguments_handle, &
    kim_compute_argument_name_partial_particle_energy, & 
    num_atoms, enepot, ierr2)
  ierr = ierr + ierr2

  if (ierr .ne. 0) then
    call kim_log_entry(model_compute_arguments_handle, & 
      kim_log_verbosity_error, "get_argument_pointer")
    return
  end if

  do cindex = 1, num_atoms
    ierr = 1 ! assumes an error
    do aux_ind = 1, buf%num_species
      if (particle_species_codes(cindex) .eq. buf%species_codes(aux_ind)) then
        ierr = 0
        exit
      end if
    end do
    if (ierr .ne. 0) then
      call kim_log_entry(model_compute_arguments_handle, & 
        kim_log_verbosity_error,  "Unexpected species code detected")
      return
    end if
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Compute energy and forces
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  !
  if (associated(energy)) then
    comp_energy =  1
  else
    comp_energy = 0
  end if
  if (associated(forces)) then
    comp_forces = 1
  else
    comp_forces = 0
  end if
  if (associated(enepot)) then
    comp_enepot = 1
  else
    comp_enepot = 0
  end if

  ierr = 0 
  allocate(neighborlists(num_atoms))
  do cindex = 1, num_atoms
    ierr = 0
    if (particle_contributing(cindex) .eq. 1) then
      call kim_get_neighbor_list(&
        model_compute_arguments_handle, 1, cindex, num_n, &
        neighbors_of_particle, ierr)
      allocate(neighborlists(cindex)%onedarray(num_n))
        
      do n_ind = 1, num_n
        neighborlists(cindex)%onedarray(n_ind) = neighbors_of_particle(n_ind)
      end do
      if (ierr .ne. 0) then
        call kim_log_entry(model_compute_arguments_handle, &
          kim_log_verbosity_error, "kim_api_get_neigh")
        ierr = 1
        return
      end if
    end if
  end do

  if (comp_enepot.eq.1) enepot = 0.0_cd
  if (comp_energy.eq.1) energy = 0.0_cd
  if (comp_forces.eq.1)  forces = 0.0_cd

  ! Initialize energy and forces
  if (comp_enepot.eq.1 .OR. comp_energy.eq.1) then
    energy = 0.0_cd
  end if
  if (comp_forces.eq.1) then
    do cindex = 1, num_atoms
      do force_dir = 1, 3
        forces(force_dir, cindex) = 0.0_cd
      end do
    end do
  end if
          
  ! Compute energy and forces directly
  if (comp_enepot.eq.1 .OR. comp_energy.eq.1 .OR. comp_forces.eq.1) then
    do cindex = 1, num_atoms
      if (particle_contributing(cindex) .eq. 1) then

        csymbol = atomic_number_codes(cindex)
        ri = coor(:, cindex)
        num_n = size(neighborlists(cindex)%onedarray)

        allocate(fingerprint(num_gs))
        allocate(ind_n(num_n))
        allocate(atomic_num_n(num_n))
        allocate(pos_n(num_n, 3))

        do n_ind = 1, num_n
          nindex = neighborlists(cindex)%onedarray(n_ind)
          ind_n(n_ind) = nindex
          atomic_num_n(n_ind) = atomic_number_codes(nindex)
          pos_n(n_ind, 1) = coor(1, nindex)
          pos_n(n_ind, 2) = coor(2, nindex)
          pos_n(n_ind, 3) = coor(3, nindex)
        end do
        
        fingerprint = calculate_fingerprint(&
          num_gs, num_blocks, buf%descriptors, &
          num_n, &
          csymbol, atomic_num_n, pos_n, &
          ri, &
          buf%cutoff(1), buf%cutofffn_code, buf%p_gamma)
  
        do ensemble_count = 1, buf%ensemble_size
          model_ind = buf%model_ensemble(ensemble_count)%navigator(csymbol)
          atom_energy = calculate_atomic_energy(&
            num_gs, &
            fingerprint, &
            buf%model_ensemble(ensemble_count)%min_fingerprints, &
            buf%model_ensemble(ensemble_count)%max_fingerprints, &
            size(buf%model_ensemble(ensemble_count)%networks(model_ind)%parameters), & 
            buf%model_ensemble(ensemble_count)%networks(model_ind)%parameters, &
            buf%model_ensemble(ensemble_count)%networks(model_ind)%no_layers, & 
            buf%model_ensemble(ensemble_count)%networks(model_ind)%no_nodes)
          energy = energy + atom_energy / buf%ensemble_size
        end do 

        ! IT IS HOW THE AFP OF THE CENTER CHANGES W.R.T ALL ITS NEIGHBORS
        ! How many neighbors the atom will affect
        ! cindex ==> n_n, nindex ==> n
        do n_ind = 1, num_n + 1
          if (n_ind .eq. 1) then
            nindex = cindex
          else
            nindex = ind_n(n_ind - 1)
          end if

          do force_dir = 0, 2
          
            allocate(fingerprintprime(num_gs))
            fingerprintprime = calculate_fingerprintprime(&
              num_gs, num_blocks, buf%descriptors, &
              num_n, &
              ind_n, atomic_num_n, pos_n, &
              buf%cutoff(1), buf%cutofffn_code, buf%p_gamma, &
              csymbol, cindex, coor(:, cindex), &
              nindex, force_dir)
              
            do ensemble_count = 1, buf%ensemble_size
              model_ind = buf%model_ensemble(ensemble_count)%navigator(csymbol)
              dforce = calculate_force(&
                num_gs, &
                fingerprint, fingerprintprime, &
                buf%model_ensemble(ensemble_count)%min_fingerprints, &
                buf%model_ensemble(ensemble_count)%max_fingerprints, &
                size(buf%model_ensemble(ensemble_count)%networks(model_ind)%parameters), &
                buf%model_ensemble(ensemble_count)%networks(model_ind)%parameters, &
                buf%model_ensemble(ensemble_count)%networks(model_ind)%no_layers, &
                buf%model_ensemble(ensemble_count)%networks(model_ind)%no_nodes)
              forces(force_dir + 1, nindex) = forces(force_dir + 1, nindex) + dforce / buf%ensemble_size
            end do 
            deallocate(fingerprintprime)
            
          end do
        end do

        deallocate(ind_n)
        deallocate(atomic_num_n)
        deallocate(pos_n)
        deallocate(fingerprint)
      end if
    end do
  end if

  do cindex = 1, num_atoms
    if (particle_contributing(cindex) .eq. 1) then
      deallocate(neighborlists(cindex)%onedarray)
    end if
  end do
  deallocate(neighborlists)

  ierr = 0
  return

end subroutine Compute_Energy_Forces

function calculate_fingerprintprime(&
  num_gs, num_blocks, gs, &
  num_n_n, &
  ind_n_n, atomic_num_n_n, pos_n_n, &
  rc, cutofffn_code, p_gamma, &
  atomic_num_n, ind_n, pos_n, & 
  ind_center, force_dir) result(derafp)
  
  implicit none

  integer(c_int):: num_gs, num_blocks, num_n_n, force_dir
  integer(c_int):: des_count, block_count
  type(atomic_descriptors), dimension(num_blocks):: gs
  
  integer(c_int):: ind_n, ind_center
  real(c_double), dimension(3):: pos_n
  
  integer(c_int):: atomic_num_n
  integer(c_int), dimension(num_n_n):: atomic_num_n_n, ind_n_n
  real(c_double), dimension(num_n_n, 3):: pos_n_n
  real(c_double)::  rc
  integer(c_int):: cutofffn_code, p_gamma

  real(c_double), dimension(num_gs):: derafp

  des_count = 1
  do block_count = 1, num_blocks
    if (gs(block_count)%symmetry_type == 'g2') then
      derafp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g2_prime(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n_n, &
          ind_n_n, atomic_num_n_n, pos_n_n, &
          ind_n, atomic_num_n, pos_n, &
          ind_center, force_dir, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes

    else if (gs(block_count)%symmetry_type == 'g3') then
      derafp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g3_prime(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n_n, &
          ind_n_n, atomic_num_n_n, pos_n_n, &
          ind_n, atomic_num_n, pos_n, &
          ind_center, force_dir, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes
      
    else if (gs(block_count)%symmetry_type == 'g3_angle') then
      derafp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g3_angle_prime(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n_n, &
          ind_n_n, atomic_num_n_n, pos_n_n, &
          ind_n, atomic_num_n, pos_n, &
          ind_center, force_dir, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%g_gamma, gs(block_count)%g_zeta, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes
      
    else if (gs(block_count)%symmetry_type == 'g3_angle_partial') then
      derafp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g3_angle_partial_prime(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n_n, &
          ind_n_n, atomic_num_n_n, pos_n_n, &
          ind_n, atomic_num_n, pos_n, &
          ind_center, force_dir, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%g_gamma, gs(block_count)%g_zeta, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes

    else if (gs(block_count)%symmetry_type == 'g4') then
      derafp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g4_prime(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n_n, &
          ind_n_n, atomic_num_n_n, pos_n_n, &
          ind_n, atomic_num_n, pos_n, &
          ind_center, force_dir, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%atomic_num_third, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes
    end if
  end do

end function calculate_fingerprintprime

function calculate_fingerprint(&
  num_gs, num_blocks, gs, &
  num_n, & 
  atomic_num, atomic_num_n, pos_n, & 
  ri, & 
  rc, cutofffn_code, p_gamma) result(afp)
  
  implicit none

  integer(c_int):: num_gs, num_blocks
  type(atomic_descriptors), dimension(num_blocks) :: gs
  integer(c_int):: des_count, block_count

  integer(c_int):: num_n, atomic_num
  integer(c_int), dimension(num_n):: atomic_num_n
  real(c_double), dimension(num_n, 3):: pos_n
  real(c_double), dimension(3):: ri
  real(c_double)::  rc
  integer(c_int):: cutofffn_code, p_gamma
  real(c_double), dimension(num_gs):: afp
          
  des_count = 1
  do block_count = 1, num_blocks
    if (gs(block_count)%symmetry_type == 'g2') then
  
      afp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g2(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n, &
          atomic_num, ri, &
          atomic_num_n, pos_n, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes

    else if (gs(block_count)%symmetry_type == 'g3') then
      afp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g3(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n, &
          atomic_num, ri, &
          atomic_num_n, pos_n, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes

    else if (gs(block_count)%symmetry_type == 'g3_angle') then
      afp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g3_angle(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n, &
          atomic_num, ri, &
          atomic_num_n, pos_n, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%g_gamma, gs(block_count)%g_zeta, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes
      
    else if (gs(block_count)%symmetry_type == 'g3_angle_partial') then
      afp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g3_angle_partial(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n, &
          atomic_num, ri, &
          atomic_num_n, pos_n, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%g_gamma, gs(block_count)%g_zeta, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes

    else if (gs(block_count)%symmetry_type == 'g4') then
      afp(des_count: des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes - 1) = &
        calculate_g4(&
          gs(block_count)%num_rscalars, gs(block_count)%r_scalars, &
          num_n, &
          atomic_num, ri, &
          atomic_num_n, pos_n, &
          gs(block_count)%atomic_num_first, &
          gs(block_count)%atomic_num_second, &
          gs(block_count)%atomic_num_third, &
          gs(block_count)%num_rtypes, &
          gs(block_count)%r_types, &
          rc, cutofffn_code, p_gamma)
      des_count = des_count + gs(block_count)%num_rscalars * gs(block_count)%num_rtypes
  
    end if
  end do

end function calculate_fingerprint

!-------------------------------------------------------------------------------
! Model driver refresh routine
!-------------------------------------------------------------------------------
recursive subroutine refresh(model_refresh_handle, ierr) bind(c)

  implicit none

  !-- Transferred variables
  type(kim_model_refresh_handle_type), intent(inout) :: model_refresh_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  real(c_double) energy_at_cutoff
  type(BUFFER_TYPE), pointer :: buf; type(c_ptr) :: pbuf

  ! get model buffer from KIM object
  call kim_get_model_buffer_pointer(model_refresh_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  call kim_set_influence_distance_pointer(model_refresh_handle, &
    buf%influence_distance(1))
  call kim_set_neighbor_list_pointers(model_refresh_handle, &
    1, buf%influence_distance, buf%not_requesting_n_from_ghost)

  ! Set new values in KIM object and buffer
  buf%influence_distance(1) = buf%cutoff(1)

  ierr = 0
  return

end subroutine refresh

!-------------------------------------------------------------------------------
! Model driver destroy routine
!-------------------------------------------------------------------------------
recursive subroutine destroy(model_destroy_handle, ierr) bind(c)
  implicit none

  !-- Transferred variables
  type(kim_model_destroy_handle_type), intent(inout) :: model_destroy_handle
  integer(c_int), intent(out) :: ierr

  !-- Local variables
  type(BUFFER_TYPE), pointer :: buf; type(c_ptr) :: pbuf

  ! get model buffer from KIM object
  call kim_get_model_buffer_pointer(model_destroy_handle, pbuf)
  call c_f_pointer(pbuf, buf)

  deallocate(buf)

  ierr = 0
  return

end subroutine destroy

!-------------------------------------------------------------------------------
! Model driver compute arguments create routine
!-------------------------------------------------------------------------------
recursive subroutine compute_arguments_create(model_compute_handle, &
  model_compute_arguments_create_handle, ierr) bind(c)
  
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_create_handle_type), intent(inout) :: &
    model_compute_arguments_create_handle
  integer(c_int), intent(out) :: ierr

  integer(c_int) :: ierr2

  ierr = 0
  ierr2 = 0

  ! register arguments
  call kim_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_energy, &
    kim_support_status_optional, ierr)
  call kim_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_forces, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  call kim_set_argument_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_argument_name_partial_particle_energy, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    call kim_log_entry(model_compute_arguments_create_handle, & 
      KIM_LOG_VERBOSITY_ERROR, & 
      "Unable to register arguments support_status")
    goto 42
  end if

  ! register callbacks
  call kim_set_callback_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_callback_name_process_dedr_term, &
    kim_support_status_optional, ierr)
  call kim_set_callback_support_status( &
    model_compute_arguments_create_handle, &
    kim_compute_callback_name_process_d2edr2_term, &
    kim_support_status_optional, ierr2)
  ierr = ierr + ierr2
  if (ierr /= 0) then
    call kim_log_entry(model_compute_arguments_create_handle, & 
      KIM_LOG_VERBOSITY_ERROR, & 
      "Unable to register callbacks support_status")
    goto 42
  end if

  ierr = 0
  42 continue
  return

end subroutine compute_arguments_create

!-------------------------------------------------------------------------------
! Model driver compute arguments destroy routine
!-------------------------------------------------------------------------------
recursive subroutine compute_arguments_destroy(model_compute_handle, &
  model_compute_arguments_destroy_handle, ierr) bind(c)
  
  implicit none

  !-- Transferred variables
  type(kim_model_compute_handle_type), intent(in) :: model_compute_handle
  type(kim_model_compute_arguments_destroy_handle_type), intent(inout) :: &
    model_compute_arguments_destroy_handle
  integer(c_int), intent(out) :: ierr

  ! nothing to be done

  ierr = 0
  return
  
end subroutine compute_arguments_destroy

end module my_driver

!-------------------------------------------------------------------------------
! Model driver create routine (REQUIRED)
!-------------------------------------------------------------------------------
recursive subroutine model_driver_create_routine(&
  model_driver_create_handle, &
  requested_length_unit, requested_energy_unit, requested_charge_unit, &
  requested_temperature_unit, requested_time_unit, ierr) bind(c)

  use, intrinsic :: iso_c_binding
  use my_driver
  use kim_model_driver_headers_module

  implicit none

  integer(c_int), parameter :: cd = c_double ! used for literal constants

  !-- Transferred variables
  type(kim_model_driver_create_handle_type), intent(inout) &
    :: model_driver_create_handle
  type(kim_length_unit_type), intent(in), value :: requested_length_unit
  type(kim_energy_unit_type), intent(in), value :: requested_energy_unit
  type(kim_charge_unit_type), intent(in), value :: requested_charge_unit
  type(kim_temperature_unit_type), intent(in), value :: & 
    requested_temperature_unit
  type(kim_time_unit_type), intent(in), value :: requested_time_unit
  integer(c_int), intent(out) :: ierr
  !-- Local variables
  integer(c_int) :: number_of_parameter_files
  character(len=1024, kind=c_char) :: parameter_file_name
  integer(c_int) :: ierr2
  type(BUFFER_TYPE), pointer :: buf
  
  integer(c_int) :: num_blocks, block_count
  integer(c_int):: num_rscalars, num_rtypes, num_descriptors
  integer(c_int) :: aux_ind, no_params

  integer(c_int) :: ensemble_count
  integer(c_int) :: aux_species, aux_navigator
  integer(c_int) :: num_models, max_atomicn, model_count
  
  ierr = 0
    
  call kim_set_units(&
    model_driver_create_handle, &
    requested_length_unit, &
    requested_energy_unit, &
    kim_charge_unit_unused, &
    kim_temperature_unit_unused, &
    kim_time_unit_unused, ierr)
  if (ierr .ne. 0) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, "Unable to set units")
    goto 42
  end if

  ! register numbering
  call kim_set_model_numbering( &
    model_driver_create_handle, kim_numbering_one_based, ierr)
  if (ierr .ne. 0) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, "Unable to set numbering")
    goto 42
  end if

  ! store callback pointers in KIM object
  call kim_set_routine_pointer( &
    model_driver_create_handle, &
    kim_model_routine_name_compute, &   
    kim_language_name_fortran, 1, &
    c_funloc(Compute_Energy_Forces), ierr)
  call kim_set_routine_pointer( &
    model_driver_create_handle, &
    kim_model_routine_name_compute_arguments_create, &
    kim_language_name_fortran, 1, &
    c_funloc(compute_arguments_create), ierr2)
  ierr = ierr + ierr2
  call kim_set_routine_pointer( &
    model_driver_create_handle, &
    kim_model_routine_name_refresh, &
    kim_language_name_fortran, 1, &
    c_funloc(refresh), ierr2)
  ierr = ierr + ierr2
  call kim_set_routine_pointer( &
    model_driver_create_handle, &
    kim_model_routine_name_compute_arguments_destroy, &
    kim_language_name_fortran, 1, &
    c_funloc(compute_arguments_destroy), ierr2)
  ierr = ierr + ierr2
  call kim_set_routine_pointer( &
    model_driver_create_handle, &
    kim_model_routine_name_destroy, &
    kim_language_name_fortran, 1, &
    c_funloc(destroy), ierr2)
  ierr = ierr + ierr2
  if (ierr .ne. 0) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, "Unable to store callback pointers")
    ierr = 1
    goto 42
  end if

  ! process parameter files
  call kim_get_number_of_parameter_files( &
    model_driver_create_handle, number_of_parameter_files)
  if (number_of_parameter_files .ne. 1) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, "Wrong number of parameter files")
    ierr = 1
    goto 42
  end if

  ! allocate model_buffer object and register it in the model_drier_create object
  allocate(buf)
  call kim_set_model_buffer_pointer( &
    model_driver_create_handle, c_loc(buf))

  ! Read in model parameters from parameter file
  call kim_get_parameter_file_name( &
    model_driver_create_handle, 1, parameter_file_name, ierr)
  if (ierr .ne. 0) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, "Unable to get parameter file name")
    ierr = 1
    goto 42
  end if
  ! 10 only indicates the file
  open(10, file=parameter_file_name, status="old")
  ! for register the atomic species later 
  read(10, *, iostat=ierr, err=100) buf%num_species
  allocate(buf%species_name_strings(buf%num_species))
  allocate(buf%atomic_numbers(buf%num_species))
  read(10, *, iostat=ierr, err=100) (buf%species_name_strings(aux_ind), aux_ind = 1, buf%num_species)
  read(10, *, iostat=ierr, err=100) (buf%atomic_numbers(aux_ind), aux_ind = 1, buf%num_species)
  max_atomicn = maxval(buf%atomic_numbers)
  
  read(10, *, iostat=ierr, err=100) buf%num_blocks
  
  num_blocks = buf%num_blocks
  allocate(buf%descriptors(num_blocks))

  num_descriptors = 0
  do block_count = 1, num_blocks
    read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%symmetry_type

    if (buf%descriptors(block_count)%symmetry_type == 'g2') then
      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%atomic_num_first
      
      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rscalars
      num_rscalars = buf%descriptors(block_count)%num_rscalars
      allocate(buf%descriptors(block_count)%r_scalars(num_rscalars))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_scalars(aux_ind), aux_ind = 1, num_rscalars)

      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rtypes
      num_rtypes = buf%descriptors(block_count)%num_rtypes
      allocate(buf%descriptors(block_count)%r_types(num_rtypes))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_types(aux_ind), aux_ind = 1, num_rtypes)

      num_descriptors = num_descriptors + num_rscalars * num_rtypes
    
    else if (buf%descriptors(block_count)%symmetry_type == 'g3') then
      read(10, *, iostat=ierr, err=100) &
        buf%descriptors(block_count)%atomic_num_first, &
        buf%descriptors(block_count)%atomic_num_second
      
      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rscalars
      num_rscalars = buf%descriptors(block_count)%num_rscalars
      allocate(buf%descriptors(block_count)%r_scalars(num_rscalars))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_scalars(aux_ind), aux_ind = 1, num_rscalars)

      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rtypes
      num_rtypes = buf%descriptors(block_count)%num_rtypes
      allocate(buf%descriptors(block_count)%r_types(num_rtypes))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_types(aux_ind), aux_ind = 1, num_rtypes)

      num_descriptors = num_descriptors + num_rscalars * num_rtypes

    else if (buf%descriptors(block_count)%symmetry_type == 'g3_angle') then
      read(10, *, iostat=ierr, err=100) &
        buf%descriptors(block_count)%atomic_num_first, &
        buf%descriptors(block_count)%atomic_num_second
      read(10, *, iostat=ierr, err=100) & 
        buf%descriptors(block_count)%g_gamma, & 
        buf%descriptors(block_count)%g_zeta
      
      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rscalars
      num_rscalars = buf%descriptors(block_count)%num_rscalars
      allocate(buf%descriptors(block_count)%r_scalars(num_rscalars))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_scalars(aux_ind), aux_ind = 1, num_rscalars)

      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rtypes
      num_rtypes = buf%descriptors(block_count)%num_rtypes
      allocate(buf%descriptors(block_count)%r_types(num_rtypes))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_types(aux_ind), aux_ind = 1, num_rtypes)

      num_descriptors = num_descriptors + num_rscalars * num_rtypes
    
    else if (buf%descriptors(block_count)%symmetry_type == 'g3_angle_partial') then
      read(10, *, iostat=ierr, err=100) &
        buf%descriptors(block_count)%atomic_num_first, &
        buf%descriptors(block_count)%atomic_num_second
      read(10, *, iostat=ierr, err=100) & 
        buf%descriptors(block_count)%g_gamma, &
        buf%descriptors(block_count)%g_zeta
        
      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rscalars
      num_rscalars = buf%descriptors(block_count)%num_rscalars
      allocate(buf%descriptors(block_count)%r_scalars(num_rscalars))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_scalars(aux_ind), aux_ind = 1, num_rscalars)

      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rtypes
      num_rtypes = buf%descriptors(block_count)%num_rtypes
      allocate(buf%descriptors(block_count)%r_types(num_rtypes))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_types(aux_ind), aux_ind = 1, num_rtypes)

      num_descriptors = num_descriptors + num_rscalars * num_rtypes

    else if (buf%descriptors(block_count)%symmetry_type == 'g4') then
      read(10, *, iostat=ierr, err=100) &
        buf%descriptors(block_count)%atomic_num_first, &
        buf%descriptors(block_count)%atomic_num_second, &
        buf%descriptors(block_count)%atomic_num_third
        
      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rscalars
      num_rscalars = buf%descriptors(block_count)%num_rscalars
      allocate(buf%descriptors(block_count)%r_scalars(num_rscalars))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_scalars(aux_ind), aux_ind = 1, num_rscalars)

      read(10, *, iostat=ierr, err=100) buf%descriptors(block_count)%num_rtypes
      num_rtypes = buf%descriptors(block_count)%num_rtypes
      allocate(buf%descriptors(block_count)%r_types(num_rtypes))
      read(10, *, iostat=ierr, err=100) &
        (buf%descriptors(block_count)%r_types(aux_ind), aux_ind = 1, num_rtypes)

      num_descriptors = num_descriptors + num_rscalars * num_rtypes
      
    end if

  end do

  buf%num_descriptors = num_descriptors

  read(10, *, iostat=ierr, err=100) buf%cutoff(1)! in A
  read(10, *, iostat=ierr, err=100) buf%cutofffn_code
  read(10, *, iostat=ierr, err=100) buf%p_gamma

  read(10, *, iostat=ierr, err=100) buf%ensemble_size
  allocate(buf%model_ensemble(buf%ensemble_size))

  do ensemble_count = 1, buf%ensemble_size

    read(10, *, iostat=ierr, err=100) num_models, max_atomicn
    allocate(buf%model_ensemble(ensemble_count)%networks(num_models))
    allocate(buf%model_ensemble(ensemble_count)%navigator(max_atomicn))
    do aux_ind = 1, buf%num_species
      read(10, *, iostat=ierr, err=100) aux_species, aux_navigator
      buf%model_ensemble(ensemble_count)%navigator(aux_species) = aux_navigator
    end do
  
    do model_count = 1, num_models
      ! get the total number of layers 
      read(10, *, iostat=ierr, err=100) & 
        buf%model_ensemble(ensemble_count)%networks(model_count)%no_layers
      allocate(buf%model_ensemble(ensemble_count)%networks(model_count)%no_nodes&
        (buf%model_ensemble(ensemble_count)%networks(model_count)%no_layers))
      ! get the number of nodes in each layer 
      read(10, *, iostat=ierr, err=100) &
        (buf%model_ensemble(ensemble_count)%networks(model_count)%no_nodes(aux_ind), aux_ind = 1, &
          buf%model_ensemble(ensemble_count)%networks(model_count)%no_layers)
      ! read the raveled model parameter of the neural network
      no_params = 0
      do aux_ind = 1, buf%model_ensemble(ensemble_count)%networks(model_count)%no_layers - 1
        no_params = no_params + &
        (buf%model_ensemble(ensemble_count)%networks(model_count)%no_nodes(aux_ind) + 1) * & 
          buf%model_ensemble(ensemble_count)%networks(model_count)%no_nodes(aux_ind + 1)
      end do
      allocate(buf%model_ensemble(ensemble_count)%networks(model_count)%parameters(no_params))
      read(10, *, iostat=ierr,err=100) & 
        (buf%model_ensemble(ensemble_count)%networks(model_count)%parameters(aux_ind), aux_ind = 1, no_params)
    end do
    
    allocate(buf%model_ensemble(ensemble_count)%min_fingerprints(num_descriptors))
    allocate(buf%model_ensemble(ensemble_count)%max_fingerprints(num_descriptors))
    read(10, *, iostat=ierr, err=100) &
      (buf%model_ensemble(ensemble_count)%min_fingerprints(aux_ind), aux_ind = 1, num_descriptors)
    read(10, *, iostat=ierr, err=100) &
      (buf%model_ensemble(ensemble_count)%max_fingerprints(aux_ind), aux_ind = 1, num_descriptors)
    
  end do 

  close(10)

  goto 200
  100 continue
  ierr = 1
  call kim_log_entry(model_driver_create_handle, & 
    kim_log_verbosity_error, "Unable to read parameters")
  goto 42
  200 continue
    
  ! register species
  allocate(buf%species_names(buf%num_species))
  allocate(buf%species_codes(buf%num_species))
  do aux_ind = 1, buf%num_species
    call kim_from_string(trim(buf%species_name_strings(aux_ind)), &
      buf%species_names(aux_ind))
      buf%species_codes(aux_ind) = aux_ind
    call kim_set_species_code(&
      model_driver_create_handle, buf%species_names(aux_ind), & 
      buf%species_codes(aux_ind), ierr)
  end do

  if (ierr .ne. 0) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, & 
      "Unable to set species_names or species_codes")
    goto 42
  end if

  buf%influence_distance(1) = buf%cutoff(1)
  ! to account for the ghost atoms
  buf%not_requesting_n_from_ghost = 1

  ! store model cutoff in KIM object
  call kim_set_influence_distance_pointer( &
    model_driver_create_handle, buf%influence_distance(1))
  call kim_set_neighbor_list_pointers( &
    model_driver_create_handle, 1, buf%influence_distance, &
    buf%not_requesting_n_from_ghost)
  ! end setup buffer

  ! store in model buffer
  call kim_set_model_buffer_pointer( &
    model_driver_create_handle, c_loc(buf))

  ! set pointers to parameters in KIM object
  call kim_set_parameter_pointer( &
    model_driver_create_handle, buf%cutoff, "cutoff", &
    "Cutoff distance of the model", ierr)
  
  if (ierr .ne. 0) then
    call kim_log_entry(model_driver_create_handle, & 
      kim_log_verbosity_error, "set_parameter")
     goto 42
  end if

  ierr = 0
  42 continue
  return

end subroutine model_driver_create_routine
