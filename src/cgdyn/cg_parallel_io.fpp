!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_parallel_io_md
!> @brief   Parallel I/O module
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_parallel_io_mod

  use cg_dynamics_str_mod
  use cg_dynvars_str_mod
  use cg_ensemble_str_mod
  use cg_boundary_str_mod
  use cg_pairlist_str_mod
  use cg_enefunc_str_mod
  use cg_domain_str_mod
  use fileio_data_mod
  use fileio_mod
  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! parameters
  integer,       parameter      :: PioFileHeaderMarker = 12345

  ! structures
  type, public :: s_pio_t0_info

    ! input files
    character(200)              :: input_topfile
    character(200)              :: input_parfile

    ! energy settings
    real(wp)                    :: energy_pairlistdist
    logical                     :: energy_table
    
    ! ensemble settings
    integer                     :: ensemble_type

    ! boundary settings
    real(wp)                    :: boundary_boxsizex
    real(wp)                    :: boundary_boxsizey
    real(wp)                    :: boundary_boxsizez
    real(wp)                    :: boundary_originx
    real(wp)                    :: boundary_originy
    real(wp)                    :: boundary_originz
    integer                     :: boundary_domain_xyz

    ! selection settings    
    character(1000),allocatable :: selection_group(:)

    ! restraint settings
    integer,        allocatable :: restraint_func(:)
    character(256), allocatable :: restraint_const(:)
    character(256), allocatable :: restraint_index(:)

  end type s_pio_t0_info

  ! global variables
  logical,              public  :: pio_restart  = .false.

  type(s_pio_t0_info),  public, save, target :: pio_t0_info

  logical                       :: pio_crdtrj_hdr = .true.
  logical                       :: pio_veltrj_hdr = .true.

  integer,          allocatable :: tmp_atom_idx(:)
  real(wp),         allocatable :: tmp_atom_crd(:,:)


  ! subroutines
  public  :: pio_write_domain_trj
  public  :: pio_write_domain_rst
  public  :: pio_read_domain_rst
  public  :: pio_read_domain_str
  public  :: pio_check_compatible
  public  :: pio_check_ranked_file
  public  :: pio_get_ranked_filename
  public  :: pio_open_file
  private :: pio_get_file_endian

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_write_domain_trj
  !> @brief        write trajectory for domain
  !! @authors      NT
  !! @param[in]    unit_no    : file unit number
  !! @param[in]    boundary   : boundary information
  !! @param[in]    domain     : domain information
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    snapshot   : number of snapshot
  !! @param[in]    crd_trj    : flag for coordinate or trajectory
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_write_domain_trj(unit_no,    &
                                  boundary,   &
                                  domain,     &
                                  trajectory, &
                                  snapshot,   &
                                  crd_trj)

    ! formal arguments
    integer,                 intent(in)    :: unit_no
    type(s_boundary),        intent(in)    :: boundary
    type(s_domain),          intent(in)    :: domain
    real(dp),                intent(in)    :: trajectory(:,:)
    integer,                 intent(in)    :: snapshot
    logical,                 intent(in)    :: crd_trj

    ! local variables
    real(wp)                 :: box(6)
    integer                  :: i, ix, natom_domain, start_i


    ! write header
    !
    if (      crd_trj .and. pio_crdtrj_hdr .or. &
        .not. crd_trj .and. pio_veltrj_hdr) then

      ! my rank number
      write(unit_no) my_world_rank
      ! # of all atoms
      write(unit_no) domain%num_atom_all
      ! # of snapshot
      write(unit_no) snapshot

      if (crd_trj) then
        pio_crdtrj_hdr = .false.
      else
        pio_veltrj_hdr = .false.
      end if

    end if

    if (.not. allocated(tmp_atom_idx)) &
      allocate(tmp_atom_idx(domain%num_cell_local * MaxAtom))

    if (.not. allocated(tmp_atom_crd)) &
      allocate(tmp_atom_crd(3,domain%num_cell_local * MaxAtom))

    ! write body
    !
    natom_domain = 1
    do i = 1, domain%num_cell_local
      start_i = (i-1)*MaxAtom
      do ix = 1, domain%num_atom(i)
        tmp_atom_idx(    natom_domain) = domain%id_l2g(ix+start_i)
        tmp_atom_crd(1:3,natom_domain) = trajectory(ix+start_i,1:3)
        natom_domain = natom_domain + 1
      end do
    end do

    natom_domain = natom_domain - 1

    box = 0.0_wp
    box(1) = boundary%box_size_x_ref
    box(3) = boundary%box_size_y_ref
    box(6) = boundary%box_size_z_ref

    ! box information
    write(unit_no) box
    ! # of domain atoms
    write(unit_no) natom_domain
    ! global atom indices
    write(unit_no) tmp_atom_idx(    1:natom_domain)
    ! atom coordinates
    write(unit_no) tmp_atom_crd (1:3,1:natom_domain)

    return

  end subroutine pio_write_domain_trj

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_write_domain_rst
  !> @brief        write restart information for domain
  !! @authors      NT
  !! @param[in]    filename    : restart filename 
  !! @param[in]    boundary    : boundary information
  !! @param[in]    domain      : domain information
  !! @param[in]    enefunc     : potential energy functions information
  !! @param[in]    dynvars     : dynvars information (optional)
  !! @param[in]    dynamics    : dynamics information (optional)
  !! @param[in]    t0_info     : t=0 information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_write_domain_rst(filename,    &
                                  boundary,    &
                                  domain,      &
                                  enefunc,     &
                                  dynvars,     &
                                  dynamics,    &
                                  t0_info)
                             
    ! formal arguments
    character(*),                   intent(in) :: filename
    type(s_boundary),               intent(in) :: boundary
    type(s_domain),                 intent(in) :: domain
    type(s_enefunc),                intent(in) :: enefunc
    type(s_dynvars),      optional, intent(in) :: dynvars
    type(s_dynamics),     optional, intent(in) :: dynamics
    type(s_pio_t0_info),  optional, intent(in) :: t0_info

    ! local variables
    integer               :: file
    integer               :: i, j, ncell, nvar, nvar1, nvar2, nvar3
    integer               :: start_i
    logical               :: tip4
    character,allocatable :: bytes(:)


    !$omp critical
    if (get_extension(filename) /= 'rsa') then
      call open_data(filename, IOFileDataWrite, file)
    else
      call open_data(filename, IOFileDataWriteAscii, file)
    end if
    !$omp end critical


    ! write constant values
    !
    call write_data_real &
         (file, 'ELECOEF', ELECOEF)


    ! write maximum value for allocatable variables
    !

    ! cg_domain_str
    call write_data_integer &
         (file, 'MaxAtom', MaxAtom)
    call write_data_integer &
         (file, 'MaxMove', MaxMove)

    ! cg_enefunc_str
    call write_data_integer &
         (file, 'MaxBond', MaxBond)
    call write_data_integer &
         (file, 'MaxAngl', MaxAngl)
    call write_data_integer &
         (file, 'MaxDihe', MaxDihe)
    call write_data_integer &
         (file, 'MaxImpr', MaxImpr)
    call write_data_integer &
         (file, 'MaxVsite2', MaxVsite2)
    call write_data_integer &
         (file, 'MaxVsite3', MaxVsite3)
    call write_data_integer &
         (file, 'MaxVsite3fd', MaxVsite3fd)
    call write_data_integer &
         (file, 'MaxVsite3fad', MaxVsite3fad)
    call write_data_integer &
         (file, 'MaxVsite3out', MaxVsite3out)
    call write_data_integer &
         (file, 'MaxVsite4fdn', MaxVsite4fdn)
    call write_data_integer &
         (file, 'MaxVsiten', MaxVsiten)
    call write_data_integer &
         (file, 'BondMove', BondMove)
    call write_data_integer &
         (file, 'AnglMove', AnglMove)
    call write_data_integer &
         (file, 'DiheMove', DiheMove)
    call write_data_integer &
         (file, 'ImprMove', ImprMove)
    call write_data_integer &
         (file, 'Vsite2Move', Vsite2Move)
    call write_data_integer &
         (file, 'Vsite3Move', Vsite3Move)
    call write_data_integer &
         (file, 'Vsite3fdMove', Vsite3fdMove)
    call write_data_integer &
         (file, 'Vsite3fadMove', Vsite3fadMove)
    call write_data_integer &
         (file, 'Vsite3outMove', Vsite3outMove)
    call write_data_integer &
         (file, 'Vsite4fdnMove', Vsite4fdnMove)
    call write_data_integer &
         (file, 'VsitenMove', VsitenMove)

    ! cg_pairlist_str
    call write_data_integer &
         (file, 'MaxNb15', MaxNb15)

    ! write boundary information
    !

    call write_data_real8 &
         (file, 'boundary:box_size_x', boundary%box_size_x)
    call write_data_real8 &
         (file, 'boundary:box_size_y', boundary%box_size_y)
    call write_data_real8 &
         (file, 'boundary:box_size_z', boundary%box_size_z)
    call write_data_real &
         (file, 'boundary:origin_x', boundary%origin_x)
    call write_data_real &
         (file, 'boundary:origin_y', boundary%origin_y)
    call write_data_real &
         (file, 'boundary:origin_z', boundary%origin_z)
    call write_data_real8 &
         (file, 'boundary:box_size_x_ref', boundary%box_size_x_ref)
    call write_data_real8 &
         (file, 'boundary:box_size_y_ref', boundary%box_size_y_ref)
    call write_data_real8 &
         (file, 'boundary:box_size_z_ref', boundary%box_size_z_ref)
    call write_data_integer_array &
         (file, 'boundary:num_domain', (/3/), boundary%num_domain(1:3))
    call write_data_integer &
         (file, 'boundary:num_cells_x', boundary%num_cells_x)
    call write_data_integer &
         (file, 'boundary:num_cells_y', boundary%num_cells_y)
    call write_data_integer &
         (file, 'boundary:num_cells_z', boundary%num_cells_z)

    ! write s_domain information
    !

    call write_data_integer &
         (file, 'domain:num_atom_all', domain%num_atom_all)
    call write_data_integer &
         (file, 'domain:num_cell_local', domain%num_cell_local)

    ! domain_str::Dynvar variables
    if (allocated(domain%num_atom)) then
      ncell = size(domain%num_atom(:))
    else
      ncell = 0
    end if

    call write_data_integer &
         (file, 'domain:ncell', ncell)

    call write_data_integer_array &
         (file, 'domain:num_atom', &
             (/ncell/), domain%num_atom   (1:ncell))
    call write_data_integer_array &
         (file, 'domain:num_atom_t0', &
             (/ncell/), domain%num_atom_t0(1:ncell))

    do i = 1, ncell

      start_i = (i-1)*MaxAtom

      ! 'MaxAtom' arrays
      nvar = domain%num_atom(i)
      call write_data_real8_array &
         (file, 'domain:coord', &
             (/nvar,3,1/), domain%coord      (start_i+1:start_i+nvar,1:3))
      call write_data_real8_array &
         (file, 'domain:velocity', &
             (/3,nvar,1/), domain%velocity   (start_i+1:start_i+nvar,1:3))
      call write_data_real_array &
         (file, 'domain:trans_vec', &
             (/3,nvar,1/), domain%trans_vec  (start_i+1:start_i+nvar,1:3))
      call write_data_real_array &
         (file, 'domain:charge', &
             (/nvar,1/),   domain%charge     (start_i+1:start_i+nvar))
      call write_data_real8_array &
         (file, 'domain:mass', &
             (/nvar,1/),   domain%mass       (start_i+1:start_i+nvar))
      call write_data_integer_array &
         (file, 'domain:id_l2g', &
             (/nvar,1/),   domain%id_l2g     (start_i+1:start_i+nvar))
      call write_data_integer_array &
         (file, 'domain:atom_cls_no', &
             (/nvar,1/),   domain%atom_cls_no(start_i+1:start_i+nvar))

    end do


    ! write s_enefunc information
    !

    ncell = domain%num_cell_local

    ! enefunc_str::Bond variables
    ! enefunc_str::Angl variables
    ! enefunc_str::Dihe variables
    ! enefunc_str::Impr variables

    call write_data_integer_array &
         (file, 'enefunc:num_bond', &
             (/ncell/), enefunc%num_bond       (1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_angl', &
             (/ncell/), enefunc%num_angl       (1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_dihedral', &
             (/ncell/), enefunc%num_dihedral   (1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_rb_dihedral', &
             (/ncell/), enefunc%num_rb_dihedral(1:ncell))
    call write_data_integer_array &
         (file, 'enefunc:num_improper', &
             (/ncell/), enefunc%num_improper   (1:ncell))
    call write_data_integer &
         (file, 'enefunc:notation_14types', enefunc%notation_14types)

    do i = 1, ncell

      nvar = enefunc%num_bond(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:bond_list', &
                 (/2,nvar,1/), enefunc%bond_list     (1:2, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:bond_force_const', &
                 (/nvar,1/), enefunc%bond_force_const(     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:bond_dist_min', &
                 (/nvar,1/), enefunc%bond_dist_min   (     1:nvar, i))
      end if

      nvar = enefunc%num_angl(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:angl_list', &
                 (/3,nvar,1/), enefunc%angl_list     (1:3, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:angl_force_const', &
                 (/nvar,1/), enefunc%angl_force_const(     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:angl_theta_min', &
                 (/nvar,1/), enefunc%angl_theta_min  (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:urey_force_const', &
                 (/nvar,1/), enefunc%urey_force_const (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:urey_rmin', &
                 (/nvar,1/), enefunc%urey_rmin        (     1:nvar, i))
      end if

      nvar = enefunc%num_dihedral(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:dihe_list', &
                 (/4,nvar,1/), enefunc%dihe_list     (1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:dihe_force_const', &
                 (/nvar,1/), enefunc%dihe_force_const(     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:dihe_periodicity', &
                 (/nvar,1/), enefunc%dihe_periodicity(     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:dihe_phase', &
                 (/nvar,1/), enefunc%dihe_phase      (     1:nvar, i))
      end if

      nvar = enefunc%num_rb_dihedral(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:rb_dihe_list', &
                 (/4,nvar,1/), enefunc%rb_dihe_list(1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:rb_dihe_c', &
                 (/6,nvar,1/), enefunc%rb_dihe_c   (1:6, 1:nvar, i))
      end if

      nvar = enefunc%num_improper(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:impr_list', &
                 (/4,nvar,1/), enefunc%impr_list     (1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:impr_force_const', &
                 (/nvar,1/), enefunc%impr_force_const(     1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:impr_periodicity', &
                 (/nvar,1/), enefunc%impr_periodicity(     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:impr_phase', &
                 (/nvar,1/), enefunc%impr_phase      (     1:nvar, i))
      end if

    end do

    ! enefunc_str::EneFuncRest variables

    if(allocated(enefunc%restraint_kind)) then
      nvar  = size(enefunc%restraint_kind(:))
      call write_data_integer &
           (file, 'enefunc:EneFuncReff_vs', nvar)
    else
      nvar=0
    endif

    if(allocated(enefunc%restraint_grouplist)) then
      nvar1 = size(enefunc%restraint_grouplist(:,1))
      call write_data_integer &
         (file, 'enefunc:EneFuncReff_vs1', nvar1)
    else
      nvar1=0
    endif

    nvar2 = max(int(nvar1/2),1)

    if (nvar > 0 .and. nvar1 > 0) then

      call write_data_integer &
         (file, 'enefunc:num_restraintfuncs', enefunc%num_restraintfuncs)

      call write_data_integer &
         (file, 'enefunc:num_restraintgroups', enefunc%num_restraintgroups)

      call write_data_integer &
         (file, 'enefunc:max_restraint_numgrps', enefunc%max_restraint_numgrps)

      call write_data_integer_array &
         (file, 'enefunc:restraint_kind', &
             (/nvar/), enefunc%restraint_kind(1:nvar))

      call write_data_integer_array &
         (file, 'enefunc:restraint_grouplist', &
             (/nvar1,nvar/), enefunc%restraint_grouplist(1:nvar1,1:nvar))

      call write_data_real_array &
         (file, 'enefunc:restraint_const', &
             (/4,nvar/), enefunc%restraint_const(1:4,1:nvar))

      call write_data_real_array &
         (file, 'enefunc:restraint_ref', &
             (/2,nvar/), enefunc%restraint_ref(1:2,1:nvar))

      call write_data_integer_array &
         (file, 'enefunc:restraint_funcgrp', &
             (/nvar/), enefunc%restraint_funcgrp(1:nvar))

      call write_data_integer_array &
         (file, 'enefunc:restraint_exponent_func', &
             (/nvar/), enefunc%restraint_exponent_func(1:nvar))

      call write_data_integer_array &
         (file, 'enefunc:restraint_exponent_dist', &
             (/nvar2,nvar/), enefunc%restraint_exponent_dist(1:nvar2,1:nvar))

      call write_data_integer_array &
         (file, 'enefunc:restraint_mode', &
             (/nvar/), enefunc%restraint_mode(1:nvar))

      call write_data_real_array &
         (file, 'enefunc:restraint_weight_dist', &
             (/nvar2,nvar/), enefunc%restraint_weight_dist(1:nvar2,1:nvar))

      call write_data_real_array &
         (file, 'enefunc:restraint_wcom1', &
             (/3,nvar1/), enefunc%restraint_wcom1(1:3,1:nvar1))

      call write_data_real_array &
         (file, 'enefunc:restraint_wcom2', &
             (/3,nvar1/), enefunc%restraint_wcom2(1:3,1:nvar1))

      call write_data_real_array &
         (file, 'enefunc:restraint_wdrt', &
             (/nvar1/), enefunc%restraint_wdrt(1:nvar2))

      call write_data_logical &
         (file, 'enefunc:restraint_posi', enefunc%restraint_posi)

      call write_data_logical &
         (file, 'enefunc:restraint_rmsd', enefunc%restraint_rmsd)

    end if

    ! enefunc_str::EneFuncRest variables

    call write_data_logical &
         (file, 'enefunc:restraint', enefunc%restraint)

    call write_data_integer_array &
         (file, 'enefunc:num_restraint', &
             (/ncell/), enefunc%num_restraint(1:ncell))

    do i = 1, ncell

      nvar = enefunc%num_restraint(i)
      if (nvar == 0) &
        cycle
      call write_data_integer_array &
           (file, 'enefunc:restraint_atom', &
               (/nvar,1/),    enefunc%restraint_atom (     1:nvar, i))
      call write_data_real_array &
           (file, 'enefunc:restraint_force', &
               (/4,nvar,1/),  enefunc%restraint_force(1:4, 1:nvar, i))
      call write_data_real_array &
           (file, 'enefunc:restraint_coord', &
               (/3,nvar,1/),  enefunc%restraint_coord(1:3, 1:nvar, i))

    end do

    ! enefunc_str::Vsite2    variables
    ! enefunc_str::Vsite3    variables
    ! enefunc_str::Vsite3fd  variables
    ! enefunc_str::Vsite3fad variables
    ! enefunc_str::Vsite3out variables
    ! enefunc_str::Vsite4fdn variables
    ! enefunc_str::Vsiten    variables

    call write_data_integer_array &
         (file, 'enefunc:num_vsite2', &
             (/ncell/), enefunc%num_vsite2(1:ncell))

    call write_data_integer_array &
         (file, 'enefunc:num_vsite3', &
             (/ncell/), enefunc%num_vsite3(1:ncell))

    call write_data_integer_array &
         (file, 'enefunc:num_vsite3fd', &
             (/ncell/), enefunc%num_vsite3fd(1:ncell))

    call write_data_integer_array &
         (file, 'enefunc:num_vsite3fad', &
             (/ncell/), enefunc%num_vsite3fad(1:ncell))

    call write_data_integer_array &
         (file, 'enefunc:num_vsite3out', &
             (/ncell/), enefunc%num_vsite3out(1:ncell))

    call write_data_integer_array &
         (file, 'enefunc:num_vsite4fdn', &
             (/ncell/), enefunc%num_vsite4fdn(1:ncell))

    call write_data_integer_array &
         (file, 'enefunc:num_vsiten', &
             (/ncell/), enefunc%num_vsiten(1:ncell))

    do i = 1, ncell

      nvar = enefunc%num_vsite2(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsite2_list', &
                 (/3,nvar,1/), enefunc%vsite2_list(1:3, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite2_a', &
                 (/nvar,1/), enefunc%vsite2_a     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsite3_list', &
                 (/4,nvar,1/), enefunc%vsite3_list(1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3_a', &
                 (/nvar,1/), enefunc%vsite3_a     (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3_b', &
                 (/nvar,1/), enefunc%vsite3_b     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3fd(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsite3fd_list', &
                 (/4,nvar,1/), enefunc%vsite3fd_list(1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3fd_a', &
                 (/nvar,1/), enefunc%vsite3fd_a     (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3fd_d', &
                 (/nvar,1/), enefunc%vsite3fd_d     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3fad(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsite3fad_list', &
                 (/4,nvar,1/), enefunc%vsite3fad_list(1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3fad_theta', &
                 (/nvar,1/), enefunc%vsite3fad_theta (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3fad_d', &
                 (/nvar,1/), enefunc%vsite3fad_d     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3out(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsite3out_list', &
                 (/4,nvar,1/), enefunc%vsite3out_list(1:4, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3out_a', &
                 (/nvar,1/), enefunc%vsite3out_a     (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3out_b', &
                 (/nvar,1/), enefunc%vsite3out_b     (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite3out_c', &
                 (/nvar,1/), enefunc%vsite3out_c     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite4fdn(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsite4fdn_list', &
                 (/5,nvar,1/), enefunc%vsite4fdn_list(1:5, 1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite4fdn_a', &
                 (/nvar,1/), enefunc%vsite4fdn_a     (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite4fdn_b', &
                 (/nvar,1/), enefunc%vsite4fdn_b     (     1:nvar, i))
        call write_data_real_array &
             (file, 'enefunc:vsite4fdn_c', &
                 (/nvar,1/), enefunc%vsite4fdn_c     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsiten(i)
      if (nvar > 0) then
        call write_data_integer_array &
             (file, 'enefunc:vsiten_list', &
                 (/10,nvar,1/), enefunc%vsiten_list(1:10, 1:nvar, i))
        call write_data_integer_array &
             (file, 'enefunc:vsiten_n', &
                 (/nvar,1/), enefunc%vsiten_n      (      1:nvar, i))
      end if

    end do

    ! enefunc_str::nonbond parameters
    
    nvar = enefunc%num_atom_cls

    call write_data_integer &
         (file, 'enefunc:num_atom_cls', nvar)

    if (nvar > 0) then
      call write_data_real_array &
           (file, 'enefunc:nb14_lj12', &
               (/nvar,nvar/), enefunc%nb14_lj12(1:nvar,1:nvar))
      call write_data_real_array &
           (file, 'enefunc:nb14_lj6', &
               (/nvar,nvar/), enefunc%nb14_lj6 (1:nvar,1:nvar))
      call write_data_real_array &
           (file, 'enefunc:nonb_lj12', &
               (/nvar,nvar/), enefunc%nonb_lj12(1:nvar,1:nvar))
      call write_data_real_array &
           (file, 'enefunc:nonb_lj6', &
               (/nvar,nvar/), enefunc%nonb_lj6 (1:nvar,1:nvar))
      call write_data_real_array &
           (file, 'enefunc:lj_coef', &
               (/   2,nvar/), lj_coef          (1:2,   1:nvar))
    end if

    ! enefunc_str:: AMBERScale
    if (allocated(enefunc%dihe_scnb)) then
      nvar = size(enefunc%dihe_scnb)
    else
      nvar = 0
    end if

    if (nvar > 0) then
      call write_data_real_array &
           (file, 'enefunc:dihe_scnb', &
               (/nvar/), enefunc%dihe_scnb(0:nvar-1))
      call write_data_real_array &
           (file, 'enefunc:dihe_scee', &
               (/nvar/), enefunc%dihe_scee(0:nvar-1))

    end if

    ! enefunc_str::other variables

    call write_data_real &
         (file, 'enefunc:table:fudge_lj', enefunc%fudge_lj)
    call write_data_real &
         (file, 'enefunc:table:fudge_qq', enefunc%fudge_qq)
    call write_data_integer &
         (file, 'enefunc:table:excl_level', enefunc%excl_level)

    ! dynvars information
    !

    if (present(dynvars)) then

      call write_data_real8 &
           (file, 'dynvars:thermostat_momentum', dynvars%thermostat_momentum)
      call write_data_real8_array &
           (file, 'dynvars:barostat_momentum', &
               (/3/), dynvars%barostat_momentum(1:3))

    end if


    ! dynamics information
    !

    if (present(dynamics)) then

      call write_data_integer &
           (file, 'dynamics:iseed', dynamics%iseed)

    end if


    ! random internal states
    !
    call random_stock_tobyte(bytes)

    call write_data_byte_array &
         (file, 'random', (/size(bytes)/), bytes)

    deallocate(bytes)


    ! t0 information
    !
    if (present(t0_info)) then

      call write_data_byte_array &
           (file, 't0_info:input_topfole', &
               (/200/), t0_info%input_topfile)
      call write_data_byte_array &
           (file, 't0_info:input_parfile', &
               (/200/), t0_info%input_parfile)
      call write_data_real &
           (file, 't0_info:energy_pairlistdist', &
               t0_info%energy_pairlistdist)
      call write_data_logical &
           (file, 't0_info:energy_table', &
               t0_info%energy_table)
      call write_data_integer &
           (file, 't0_info:ensemble_type', &
               t0_info%ensemble_type)
      call write_data_real &
           (file, 't0_info:boundary_boxsizex', &
               t0_info%boundary_boxsizex)
      call write_data_real &
           (file, 't0_info:boundary_boxsizey', &
               t0_info%boundary_boxsizey)
      call write_data_real &
           (file, 't0_info:boundary_boxsizez', &
               t0_info%boundary_boxsizez)
      call write_data_real &
           (file, 't0_info:boundary_originx', &
               t0_info%boundary_originx)
      call write_data_real &
           (file, 't0_info:boundary_originy', &
               t0_info%boundary_originy)
      call write_data_real &
           (file, 't0_info:boundary_originz', &
               t0_info%boundary_originz)
      call write_data_integer & 
           (file, 't0_info:boundary_domain_xyz', &
               t0_info%boundary_domain_xyz)

      if (allocated(t0_info%selection_group)) then
        nvar = size(t0_info%selection_group)
      else
        nvar = 0
      end if

      if (nvar > 0) then
        call write_data_byte_array &
             (file, 't0_info:selection_group', &
                 (/1000*nvar/), t0_info%selection_group(1:nvar))
      end if

      if (allocated(t0_info%restraint_func)) then
        nvar = size(t0_info%restraint_func)
      else
        nvar = 0
      end if

      if (nvar > 0) then
        call write_data_integer_array &
             (file, 't0_info:restraint_func', &
                 (/nvar/), t0_info%restraint_func (1:nvar))
        call write_data_byte_array &
             (file, 't0_info:restraint_const', &
                 (/nvar*256/), t0_info%restraint_const(1:nvar))
        call write_data_byte_array &
             (file, 't0_info:restraint_index', &
                 (/nvar*256/), t0_info%restraint_index(1:nvar))
      end if

    else

      call write_data_byte_array &
           (file, 't0_info:input_topfole', &
               (/200/), pio_t0_info%input_topfile)
      call write_data_byte_array &
           (file, 't0_info:input_parfile', &
               (/200/), pio_t0_info%input_parfile)
      call write_data_real &
           (file, 't0_info:energy_pairlistdist', &
               pio_t0_info%energy_pairlistdist)
      call write_data_logical &
           (file, 't0_info:energy_table', &
               pio_t0_info%energy_table)
      call write_data_integer &
           (file, 't0_info:ensemble_type', &
               pio_t0_info%ensemble_type)
      call write_data_real &
           (file, 't0_info:boundary_boxsizex', &
               pio_t0_info%boundary_boxsizex)
      call write_data_real &
           (file, 't0_info:boundary_boxsizey', &
               pio_t0_info%boundary_boxsizey)
      call write_data_real &
           (file, 't0_info:boundary_boxsizez', &
               pio_t0_info%boundary_boxsizez)
      call write_data_real &
           (file, 't0_info:boundary_originx', &
               pio_t0_info%boundary_originx)
      call write_data_real &
           (file, 't0_info:boundary_originy', &
               pio_t0_info%boundary_originy)
      call write_data_real &
           (file, 't0_info:boundary_originz', &
               pio_t0_info%boundary_originz)
      call write_data_integer & 
           (file, 't0_info:boundary_domain_xyz', &
               pio_t0_info%boundary_domain_xyz)

      if (allocated(pio_t0_info%selection_group)) then
        nvar = size(pio_t0_info%selection_group)
      else
        nvar = 0
      end if

      if (nvar > 0) then
        call write_data_byte_array &
             (file, 't0_info:selection_group', &
                 (/1000*nvar/), pio_t0_info%selection_group(1:nvar))
      end if

      if (allocated(pio_t0_info%restraint_func)) then
        nvar = size(pio_t0_info%restraint_func)
      else
        nvar = 0
      end if

      if (nvar > 0) then

        call write_data_integer_array &
             (file, 't0_info:restraint_func', &
                 (/nvar/), pio_t0_info%restraint_func (1:nvar))
        call write_data_byte_array &
             (file, 't0_info:restraint_const', &
                 (/nvar*256/), pio_t0_info%restraint_const(1:nvar))
        call write_data_byte_array &
             (file, 't0_info:restraint_index', &
                 (/nvar*256/), pio_t0_info%restraint_index(1:nvar))
      end if

    end if


    ! write restart information
    !       (This domain restart file is RESTART or not.)
    call write_data_logical &
         (file, 'pio_restart', pio_restart)


    !$omp critical
    call close_data(file)
    !$omp end critical

    return

  end subroutine pio_write_domain_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_read_domain_rst
  !> @brief        read restart information for domain
  !! @authors      NT
  !! @param[in]    filename    : restart filename 
  !! @param[inout] boundary    : boundary information
  !! @param[inout] domain      : domain information
  !! @param[inout] enefunc     : potential energy functions information
  !! @param[inout] dynvars     : dynvars information (optional)
  !! @param[inout] dynamics    : dynamics information (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_domain_rst(filename,    &
                                 boundary,    &
                                 domain,      &
                                 enefunc,     &
                                 dynvars,     &
                                 dynamics)
                             
    ! formal arguments
    character(*),                   intent(in)    :: filename
    type(s_boundary),               intent(inout) :: boundary
    type(s_domain),                 intent(inout) :: domain
    type(s_enefunc),                intent(inout) :: enefunc
    type(s_dynvars),      optional, intent(inout) :: dynvars
    type(s_dynamics),     optional, intent(inout) :: dynamics

    ! local variables
    integer               :: file
    integer               :: i, j, ncell, nvar, nvar1, nvar2, nvar3
    integer               :: start_i, start_icel
    integer               :: ix, icel, ig, nbyte
    integer               :: iseed_tmp = 0
    logical               :: tip4
    character,allocatable :: bytes(:)


    ! initialize
    !
    call init_boundary(boundary)
    call init_domain  (domain)
    call init_enefunc (enefunc)
    if (present(dynvars))     call init_dynvars    (dynvars)
    if (present(dynamics))    call init_dynamics   (dynamics)


    call open_data(filename, IOFileDataRead, file)


    ! read constant values
    !
    call read_data_real &
         (file, 'ELECOEF', ELECOEF)

    ! read maximum value for allocatable variables
    !

    ! cg_domain_str
    call read_data_integer &
         (file, 'MaxAtom', MaxAtom)
    call read_data_integer &
         (file, 'MaxMove', MaxMove)

    ! cg_enefunc_str
    call read_data_integer &
         (file, 'MaxBond', MaxBond)
    call read_data_integer &
         (file, 'MaxAngl', MaxAngl)
    call read_data_integer &
         (file, 'MaxDihe', MaxDihe)
    call read_data_integer &
         (file, 'MaxImpr', MaxImpr)
    call read_data_integer &
         (file, 'MaxVsite2', MaxVsite2)
    call read_data_integer &
         (file, 'MaxVsite3', MaxVsite3)
    call read_data_integer &
         (file, 'MaxVsite3fd', MaxVsite3fd)
    call read_data_integer &
         (file, 'MaxVsite3fad', MaxVsite3fad)
    call read_data_integer &
         (file, 'MaxVsite3out', MaxVsite3out)
    call read_data_integer &
         (file, 'MaxVsite4fdn', MaxVsite4fdn)
    call read_data_integer &
         (file, 'MaxVsiten', MaxVsiten)
    call read_data_integer &
         (file, 'BondMove', BondMove)
    call read_data_integer &
         (file, 'AnglMove', AnglMove)
    call read_data_integer &
         (file, 'DiheMove', DiheMove)
    call read_data_integer &
         (file, 'ImprMove', ImprMove)
    call read_data_integer &
         (file, 'Vsite2Move', Vsite2Move)
    call read_data_integer &
         (file, 'Vsite3Move', Vsite3Move)
    call read_data_integer &
         (file, 'Vsite3fdMove', Vsite3fdMove)
    call read_data_integer &
         (file, 'Vsite3fadMove', Vsite3fadMove)
    call read_data_integer &
         (file, 'Vsite3outMove', Vsite3outMove)
    call read_data_integer &
         (file, 'Vsite4fdnMove', Vsite4fdnMove)
    call read_data_integer &
         (file, 'VsitenMove', VsitenMove)

    ! cg_pairlist_str
    call read_data_integer &
         (file, 'MaxNb15', MaxNb15)

    ! read boundary information
    !

    call read_data_real8 &
         (file, 'boundary:box_size_x', boundary%box_size_x)
    call read_data_real8 &
         (file, 'boundary:box_size_y', boundary%box_size_y)
    call read_data_real8 &
         (file, 'boundary:box_size_z', boundary%box_size_z)
    call read_data_real &
         (file, 'boundary:origin_x', boundary%origin_x)
    call read_data_real &
         (file, 'boundary:origin_y', boundary%origin_y)
    call read_data_real &
         (file, 'boundary:origin_z', boundary%origin_z)
    call read_data_real8 &
         (file, 'boundary:box_size_x_ref', boundary%box_size_x_ref)
    call read_data_real8 &
         (file, 'boundary:box_size_y_ref', boundary%box_size_y_ref)
    call read_data_real8 &
         (file, 'boundary:box_size_z_ref', boundary%box_size_z_ref)
    call read_data_integer_array &
         (file, 'boundary:num_domain', (/3/), boundary%num_domain(1:3))
    call read_data_integer &
         (file, 'boundary:num_cells_x', boundary%num_cells_x)
    call read_data_integer &
         (file, 'boundary:num_cells_y', boundary%num_cells_y)
    call read_data_integer &
         (file, 'boundary:num_cells_z', boundary%num_cells_z)

    ! read s_domain information
    !

    call read_data_integer &
         (file, 'domain:num_atom_all', domain%num_atom_all)
    call read_data_integer &
         (file, 'domain:num_cell_local', domain%num_cell_local)

    ! domain_str::Dynvar variables
    call read_data_integer &
         (file, 'domain:ncell', ncell)
    call alloc_domain(domain, DomainDynvar, ncell, 1, 1)

    call read_data_integer_array &
         (file, 'domain:num_atom', &
             (/ncell/), domain%num_atom   (1:ncell))
    call read_data_integer_array &
         (file, 'domain:num_atom_t0', &
             (/ncell/), domain%num_atom_t0(1:ncell))

    call alloc_domain(domain, DomainDynvar_Atom, ncell, nvar1, 1)

    do i = 1, ncell

      start_i = (i-1)*MaxAtom
      ! 'MaxAtom' arrays
      nvar = domain%num_atom(i)
      call read_data_real8_array &
         (file, 'domain:coord', &
             (/nvar,3,1/), domain%coord      (start_i+1:start_i+nvar,1:3))
      call read_data_real8_array &
         (file, 'domain:velocity', &
             (/nvar,3,1/), domain%velocity   (start_i+1:start_i+nvar,1:3))
      call read_data_real_array &
         (file, 'domain:trans_vec', &
             (/nvar,3,1/), domain%trans_vec  (start_i+1:start_i+nvar,1:3))
      call read_data_real_array &
         (file, 'domain:charge', &
             (/nvar,1/),   domain%charge     (start_i+1:start_i+nvar))
      call read_data_real8_array &
         (file, 'domain:mass', &
             (/nvar,1/),   domain%mass       (start_i+1:start_i+nvar))
      call read_data_integer_array &
         (file, 'domain:id_l2g', &
             (/nvar,1/),   domain%id_l2g     (start_i+1:start_i+nvar))
      call read_data_integer_array &
         (file, 'domain:atom_cls_no', &
             (/nvar,1/),   domain%atom_cls_no(start_i+1:start_i+nvar))

    end do

    ! domain_str::Global variables
    call alloc_domain(domain, DomainGlobal, domain%num_atom_all, 1, 1)
    !! Re-compute id_g2l by id_l2g
    do icel = 1, ncell
      start_icel = (icel-1)*MaxAtom
      do ix = 1, domain%num_atom(icel)
        ig = domain%id_l2g(ix+start_icel)
        domain%id_g2l(ig) = ix+icel*MaxAtom
      end do
    end do
    

    ! read s_enefunc information
    !

    ncell = domain%num_cell_local

    ! enefunc_str::Bond variables
    ! enefunc_str::Angl variables
    ! enefunc_str::Dihe variables
    ! enefunc_str::RBDihe variables
    ! enefunc_str::Impr variables

    call alloc_enefunc(enefunc, EneFuncBase,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncBond,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncAngl,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncDihe,   ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncRBDihe, ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncImpr,   ncell, ncell)

    call read_data_integer_array &
         (file, 'enefunc:num_bond', &
             (/ncell/), enefunc%num_bond       (1:ncell))
    call read_data_integer_array &
         (file, 'enefunc:num_angl', &
             (/ncell/), enefunc%num_angl      (1:ncell))
    call read_data_integer_array &
         (file, 'enefunc:num_dihedral', &
             (/ncell/), enefunc%num_dihedral   (1:ncell))
    call read_data_integer_array &
         (file, 'enefunc:num_rb_dihedral', &
             (/ncell/), enefunc%num_rb_dihedral(1:ncell))
    call read_data_integer_array &
         (file, 'enefunc:num_improper', &
             (/ncell/), enefunc%num_improper   (1:ncell))
    call read_data_integer &
         (file, 'enefunc:notation_14types', enefunc%notation_14types)

    do i = 1, ncell

      nvar = enefunc%num_bond(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:bond_list', &
                 (/2,nvar,1/), enefunc%bond_list     (1:2, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:bond_force_const', &
                 (/nvar,1/), enefunc%bond_force_const(     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:bond_dist_min', &
                 (/nvar,1/), enefunc%bond_dist_min   (     1:nvar, i))
      end if

      nvar = enefunc%num_angl(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:angl_list', &
                 (/3,nvar,1/), enefunc%angl_list     (1:3, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:angl_force_const', &
                 (/nvar,1/), enefunc%angl_force_const(     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:angl_theta_min', &
                 (/nvar,1/), enefunc%angl_theta_min  (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:urey_force_const', &
                 (/nvar,1/), enefunc%urey_force_const (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:urey_rmin', &
                 (/nvar,1/), enefunc%urey_rmin        (     1:nvar, i))
      end if

      nvar = enefunc%num_dihedral(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:dihe_list', &
                 (/4,nvar,1/), enefunc%dihe_list     (1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:dihe_force_const', &
                 (/nvar,1/), enefunc%dihe_force_const(     1:nvar, i))
        call read_data_integer_array &
             (file, 'enefunc:dihe_periodicity', &
                 (/nvar,1/), enefunc%dihe_periodicity(     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:dihe_phase', &
                 (/nvar,1/), enefunc%dihe_phase      (     1:nvar, i))
      end if

      nvar = enefunc%num_rb_dihedral(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:rb_dihe_list', &
                 (/4,nvar,1/), enefunc%rb_dihe_list(1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:rb_dihe_c', &
                 (/6,nvar,1/), enefunc%rb_dihe_c   (1:6, 1:nvar, i))
      end if

      nvar = enefunc%num_improper(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:impr_list', &
                 (/4,nvar,1/), enefunc%impr_list     (1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:impr_force_const', &
                 (/nvar,1/), enefunc%impr_force_const(     1:nvar, i))
        call read_data_integer_array &
             (file, 'enefunc:impr_periodicity', &
                 (/nvar,1/), enefunc%impr_periodicity(     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:impr_phase', &
                 (/nvar,1/), enefunc%impr_phase      (     1:nvar, i))
      end if

    end do

    ! enefunc_str::EneFuncReff variables

    nvar = 0
    call read_data_integer &
         (file, 'enefunc:EneFuncReff_vs', nvar, .false.)

    nvar1 = 0
    call read_data_integer &
         (file, 'enefunc:EneFuncReff_vs1', nvar1, .false.)

    nvar2 = max(int(nvar1/2),1)

    if (nvar > 0 .and. nvar1 > 0) then

      call read_data_integer &
         (file, 'enefunc:num_restraintfuncs', enefunc%num_restraintfuncs)

      call read_data_integer &
         (file, 'enefunc:num_restraintgroups', enefunc%num_restraintgroups)

      call read_data_integer &
         (file, 'enefunc:max_restraint_numgrps', enefunc%max_restraint_numgrps)

      call alloc_enefunc(enefunc, EneFuncReff, nvar, nvar1)

      call read_data_integer_array &
         (file, 'enefunc:restraint_kind', &
             (/nvar/), enefunc%restraint_kind(1:nvar))

      call read_data_integer_array &
         (file, 'enefunc:restraint_grouplist', &
             (/nvar1,nvar/), enefunc%restraint_grouplist(1:nvar1,1:nvar))

      call read_data_real_array &
         (file, 'enefunc:restraint_const', &
             (/4,nvar/), enefunc%restraint_const(1:4,1:nvar))

      call read_data_real_array &
         (file, 'enefunc:restraint_ref', &
             (/2,nvar/), enefunc%restraint_ref(1:2,1:nvar))

      call read_data_integer_array &
         (file, 'enefunc:restraint_funcgrp', &
             (/nvar/), enefunc%restraint_funcgrp(1:nvar))

      call read_data_integer_array &
         (file, 'enefunc:restraint_exponent_func', &
             (/nvar/), enefunc%restraint_exponent_func(1:nvar))

      call read_data_integer_array &
         (file, 'enefunc:restraint_exponent_dist', &
             (/nvar2,nvar/), enefunc%restraint_exponent_dist(1:nvar2,1:nvar))

      call read_data_integer_array &
         (file, 'enefunc:restraint_mode', &
             (/nvar/), enefunc%restraint_mode(1:nvar))

      call read_data_real_array &
         (file, 'enefunc:restraint_weight_dist', &
             (/nvar2,nvar/), enefunc%restraint_weight_dist(1:nvar2,1:nvar))

      call read_data_real_array &
         (file, 'enefunc:restraint_wcom1', &
             (/3,nvar1/), enefunc%restraint_wcom1(1:3,1:nvar1))

      call read_data_real_array &
         (file, 'enefunc:restraint_wcom2', &
             (/3,nvar1/), enefunc%restraint_wcom2(1:3,1:nvar1))

      call read_data_real_array &
         (file, 'enefunc:restraint_wdrt', &
             (/nvar1/), enefunc%restraint_wdrt(1:nvar2))

      call read_data_logical &
         (file, 'enefunc:restraint_posi', enefunc%restraint_posi)

      call read_data_logical &
         (file, 'enefunc:restraint_rmsd', enefunc%restraint_rmsd)

    end if

    ! enefunc_str::EneFuncRest variables

    call read_data_logical &
         (file, 'enefunc:restraint', enefunc%restraint)

    call read_data_integer_array &
         (file, 'enefunc:num_restraint', &
             (/ncell/), enefunc%num_restraint(1:ncell))

    call alloc_enefunc(enefunc, EneFuncRest, ncell, ncell)

    do i = 1, ncell

      nvar = enefunc%num_restraint(i)
      if (nvar == 0) &
        cycle
      call read_data_integer_array &
           (file, 'enefunc:restraint_atom', &
               (/nvar,1/),    enefunc%restraint_atom (     1:nvar, i))
      call read_data_real_array &
           (file, 'enefunc:restraint_force', &
               (/4,nvar,1/),  enefunc%restraint_force(1:4, 1:nvar, i))
      call read_data_real_array &
           (file, 'enefunc:restraint_coord', &
               (/3,nvar,1/),  enefunc%restraint_coord(1:3, 1:nvar, i))

    end do

    ! enefunc_str::Vsite2    variables
    ! enefunc_str::Vsite3    variables
    ! enefunc_str::Vsite3fd  variables
    ! enefunc_str::Vsite3fad variables
    ! enefunc_str::Vsite3out variables
    ! enefunc_str::Vsite4fdn variables
    ! enefunc_str::Vsiten    variables

    call alloc_enefunc(enefunc, EneFuncVsite2,    ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncVsite3,    ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncVsite3fd,  ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncVsite3fad, ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncVsite3out, ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncVsite4fdn, ncell, ncell)
    call alloc_enefunc(enefunc, EneFuncVsiten,    ncell, ncell)

    call read_data_integer_array &
         (file, 'enefunc:num_vsite2', &
             (/ncell/), enefunc%num_vsite2(1:ncell), .false.)

    call read_data_integer_array &
         (file, 'enefunc:num_vsite3', &
             (/ncell/), enefunc%num_vsite3(1:ncell), .false.)

    call read_data_integer_array &
         (file, 'enefunc:num_vsite3fd', &
             (/ncell/), enefunc%num_vsite3fd(1:ncell), .false.)

    call read_data_integer_array &
         (file, 'enefunc:num_vsite3fad', &
             (/ncell/), enefunc%num_vsite3fad(1:ncell), .false.)

    call read_data_integer_array &
         (file, 'enefunc:num_vsite3out', &
             (/ncell/), enefunc%num_vsite3out(1:ncell), .false.)

    call read_data_integer_array &
         (file, 'enefunc:num_vsite4fdn', &
             (/ncell/), enefunc%num_vsite4fdn(1:ncell), .false.)

    call read_data_integer_array &
         (file, 'enefunc:num_vsiten', &
             (/ncell/), enefunc%num_vsiten(1:ncell), .false.)

    do i = 1, ncell

      nvar = enefunc%num_vsite2(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsite2_list', &
                 (/3,nvar,1/), enefunc%vsite2_list(1:3, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite2_a', &
                 (/nvar,1/), enefunc%vsite2_a     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsite3_list', &
                 (/4,nvar,1/), enefunc%vsite3_list(1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3_a', &
                 (/nvar,1/), enefunc%vsite3_a     (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3_b', &
                 (/nvar,1/), enefunc%vsite3_b     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3fd(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsite3fd_list', &
                 (/4,nvar,1/), enefunc%vsite3fd_list(1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3fd_a', &
                 (/nvar,1/), enefunc%vsite3fd_a     (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3fd_d', &
                 (/nvar,1/), enefunc%vsite3fd_d     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3fad(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsite3fad_list', &
                 (/4,nvar,1/), enefunc%vsite3fad_list(1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3fad_theta', &
                 (/nvar,1/), enefunc%vsite3fad_theta (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3fad_d', &
                 (/nvar,1/), enefunc%vsite3fad_d     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite3out(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsite3out_list', &
                 (/4,nvar,1/), enefunc%vsite3out_list(1:4, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3out_a', &
                 (/nvar,1/), enefunc%vsite3out_a     (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3out_b', &
                 (/nvar,1/), enefunc%vsite3out_b     (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite3out_c', &
                 (/nvar,1/), enefunc%vsite3out_c     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsite4fdn(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsite4fdn_list', &
                 (/5,nvar,1/), enefunc%vsite4fdn_list(1:5, 1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite4fdn_a', &
                 (/nvar,1/), enefunc%vsite4fdn_a     (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite4fdn_b', &
                 (/nvar,1/), enefunc%vsite4fdn_b     (     1:nvar, i))
        call read_data_real_array &
             (file, 'enefunc:vsite4fdn_c', &
                 (/nvar,1/), enefunc%vsite4fdn_c     (     1:nvar, i))
      end if

      nvar = enefunc%num_vsiten(i)
      if (nvar > 0) then
        call read_data_integer_array &
             (file, 'enefunc:vsiten_list', &
                 (/10,nvar,1/), enefunc%vsiten_list(1:10, 1:nvar, i))
        call read_data_integer_array &
             (file, 'enefunc:vsiten_n', &
                 (/nvar,1/), enefunc%vsiten_n      (      1:nvar, i))
      end if

    end do

    ! enefunc_str::nonbond parameters
    
    call read_data_integer &
         (file, 'enefunc:num_atom_cls', enefunc%num_atom_cls)

    nvar = enefunc%num_atom_cls

    call alloc_enefunc(enefunc, EneFuncNbon, nvar)

    if (nvar > 0) then
      call read_data_real_array &
           (file, 'enefunc:nb14_lj12', &
               (/nvar,nvar/), enefunc%nb14_lj12(1:nvar,1:nvar))
      call read_data_real_array &
           (file, 'enefunc:nb14_lj6', &
               (/nvar,nvar/), enefunc%nb14_lj6 (1:nvar,1:nvar))
      call read_data_real_array &
           (file, 'enefunc:nonb_lj12', &
               (/nvar,nvar/), enefunc%nonb_lj12(1:nvar,1:nvar))
      call read_data_real_array &
           (file, 'enefunc:nonb_lj6', &
               (/nvar,nvar/), enefunc%nonb_lj6 (1:nvar,1:nvar))
      call read_data_real_array &
           (file, 'enefunc:lj_coef', &
               (/   2,nvar/), lj_coef          (1:2,   1:nvar))
    end if

    ! enefunc_str:: AMBERScale
    call get_data_size &
         (file, 'enefunc:dihe_scnb', nvar, .false.)

    call alloc_enefunc(enefunc, EneFuncAMBERScale, nvar)

    if (nvar > 0) then
      call read_data_real_array &
           (file, 'enefunc:dihe_scnb', &
               (/nvar/), enefunc%dihe_scnb(0:nvar-1))
      call read_data_real_array &
           (file, 'enefunc:dihe_scee', &
               (/nvar/), enefunc%dihe_scee(0:nvar-1))
    end if

    ! enefunc_str::other variables

    call read_data_real &
         (file, 'enefunc:table:fudge_lj', enefunc%fudge_lj)
    call read_data_real &
         (file, 'enefunc:table:fudge_qq', enefunc%fudge_qq)
    call read_data_integer &
         (file, 'enefunc:table:excl_level', enefunc%excl_level)

    ! dynvars
    !

    if (present(dynvars)) then

      call read_data_real8 &
           (file, 'dynvars:thermostat_momentum', dynvars%thermostat_momentum)
      call read_data_real8_array &
           (file, 'dynvars:barostat_momentum', &
               (/3/), dynvars%barostat_momentum(1:3))
    end if


    ! dynamics
    !

    if (present(dynamics)) then

      call read_data_integer &
           (file, 'dynamics:iseed', dynamics%iseed)
      iseed_tmp = dynamics%iseed

    end if


    ! random internal states
    !
    call get_data_size &
         (file, 'random', nbyte, .false.)
    allocate(bytes(nbyte))

    call read_data_byte_array &
         (file, 'random', (/nbyte/), bytes(1:nbyte))

    call random_init(iseed_tmp)
    call random_stock_frombyte(bytes)
    call random_pull_stock
    
    deallocate(bytes)


    ! t0 information
    !

    call read_data_byte_array &
         (file, 't0_info:input_topfole', &
             (/200/), pio_t0_info%input_topfile)
    call read_data_byte_array &
         (file, 't0_info:input_parfile', &
             (/200/), pio_t0_info%input_parfile)
    call read_data_real &
         (file, 't0_info:energy_pairlistdist', &
             pio_t0_info%energy_pairlistdist)
    call read_data_logical &
         (file, 't0_info:energy_table', &
             pio_t0_info%energy_table)
    call read_data_integer &
         (file, 't0_info:ensemble_type', &
             pio_t0_info%ensemble_type)
    call read_data_real &
         (file, 't0_info:boundary_boxsizex', &
             pio_t0_info%boundary_boxsizex)
    call read_data_real &
         (file, 't0_info:boundary_boxsizey', &
             pio_t0_info%boundary_boxsizey)
    call read_data_real &
         (file, 't0_info:boundary_boxsizez', &
             pio_t0_info%boundary_boxsizez)
    call read_data_real &
         (file, 't0_info:boundary_originx', &
             pio_t0_info%boundary_originx)
    call read_data_real &
         (file, 't0_info:boundary_originy', &
             pio_t0_info%boundary_originy)
    call read_data_real &
         (file, 't0_info:boundary_originz', &
             pio_t0_info%boundary_originz)
    call read_data_integer & 
         (file, 't0_info:boundary_domain_xyz', &
             pio_t0_info%boundary_domain_xyz)

    call get_data_size &
         (file, 't0_info:selection_group', nvar, .false.)
    nvar = nvar / 1000
    if (allocated(pio_t0_info%selection_group)) &
      deallocate(pio_t0_info%selection_group)
    allocate(pio_t0_info%selection_group(nvar))

    if (nvar > 0) then
      call read_data_byte_array &
           (file, 't0_info:selection_group', &
               (/1000*nvar/), pio_t0_info%selection_group(1:nvar))
    end if

    call get_data_size &
         (file, 't0_info:restraint_func', nvar, .false.)
    if (allocated(pio_t0_info%restraint_func))  &
      deallocate(pio_t0_info%restraint_func,    &
                 pio_t0_info%restraint_const,   &
                 pio_t0_info%restraint_index)
    allocate(pio_t0_info%restraint_func(nvar),  &
             pio_t0_info%restraint_const(nvar), &
             pio_t0_info%restraint_index(nvar))

    if (nvar > 0) then
      call read_data_integer_array &
           (file, 't0_info:restraint_func', &
               (/nvar/), pio_t0_info%restraint_func (1:nvar))
      call read_data_byte_array &
           (file, 't0_info:restraint_const', &
               (/nvar*256/), pio_t0_info%restraint_const(1:nvar))
      call read_data_byte_array &
           (file, 't0_info:restraint_index', &
               (/nvar*256/), pio_t0_info%restraint_index(1:nvar))
    end if


    ! read restart information
    !       (Domain restart file is START or RESTART)
    call read_data_logical &
         (file, 'pio_restart', pio_restart)


    call close_data(file)

    return

  end subroutine pio_read_domain_rst

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_read_domain_str
  !> @brief        read domain information
  !! @authors      NT
  !! @param[in]    filename : restart filename 
  !! @param[inout] boundary : boundary information
  !! @param[inout] domain   : domain information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_read_domain_str(filename, boundary, domain, convert)
                             
    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_boundary),        intent(inout) :: boundary
    type(s_domain),          intent(inout) :: domain
    logical,                 intent(in)    :: convert

    ! local variables
    integer                  :: file
    integer                  :: i, ncell, nvar
    integer                  :: ix, icel, ig
    integer                  :: start_i, start_icel


    call open_data(filename, IOFileDataRead, file)


    ! read constant values
    !
    if (.not. convert) &
      call read_data_real &
         (file, 'ELECOEF', ELECOEF)


    ! read boundary information
    !

    call read_data_real8 &
         (file, 'boundary:box_size_x', boundary%box_size_x)
    call read_data_real8 &
         (file, 'boundary:box_size_y', boundary%box_size_y)
    call read_data_real8 &
         (file, 'boundary:box_size_z', boundary%box_size_z)
    call read_data_real &
         (file, 'boundary:origin_x', boundary%origin_x)
    call read_data_real &
         (file, 'boundary:origin_y', boundary%origin_y)
    call read_data_real &
         (file, 'boundary:origin_z', boundary%origin_z)
    call read_data_real8 &
         (file, 'boundary:box_size_x_ref', boundary%box_size_x_ref)
    call read_data_real8 &
         (file, 'boundary:box_size_y_ref', boundary%box_size_y_ref)
    call read_data_real8 &
         (file, 'boundary:box_size_z_ref', boundary%box_size_z_ref)
    call read_data_integer_array &
         (file, 'boundary:num_domain', (/3/), boundary%num_domain(1:3))
    call read_data_integer &
         (file, 'boundary:num_cells_x', boundary%num_cells_x)
    call read_data_integer &
         (file, 'boundary:num_cells_y', boundary%num_cells_y)
    call read_data_integer &
         (file, 'boundary:num_cells_z', boundary%num_cells_z)


    ! read s_domain information
    !

    call read_data_integer &
         (file, 'domain:num_atom_all', domain%num_atom_all)
    call read_data_integer &
         (file, 'domain:num_cell_local', domain%num_cell_local)

    ! domain_str::Dynvar variables
    call read_data_integer &
         (file, 'domain:ncell', ncell)
    call alloc_domain(domain, DomainDynvar, ncell, 1, 1)

    call read_data_integer_array &
         (file, 'domain:num_atom', &
             (/ncell/), domain%num_atom   (1:ncell))

    do i = 1, ncell

      ! 'MaxAtom' arrays
      nvar = domain%num_atom(i)
      call read_data_real8_array &
         (file, 'domain:coord', &
             (/nvar,3,1/), domain%coord   (start_i+1:start_i+nvar,1:3))
      call read_data_real8_array &
         (file, 'domain:velocity', &
             (/nvar,3,1/), domain%velocity(start_i+1:start_i+nvar,1:3))
      call read_data_integer_array &
         (file, 'domain:id_l2g', &
             (/nvar,1/),   domain%id_l2g  (start_i+1:start_i+nvar))

    end do

    call close_data(file)

    return

  end subroutine pio_read_domain_str

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_check_compatible
  !> @brief        check md input parameters are compatible between t=0 and t=n
  !! @authors      NT
  !! @param[in]    energy_pairlistdist   :      [ENERGY] pairlistdist
  !! @param[in]    energy_table          :      [ENERGY] table
  !! @param[in]    ensemble_type         :    [ENSEMBLE] ensemble
  !! @param[in]    boundary_boxsizex     :    [BOUNDARY] boxsizex
  !! @param[in]    boundary_boxsizey     :    [BOUNDARY] boxsizey
  !! @param[in]    boundary_boxsizez     :    [BOUNDARY] boxsizez
  !! @param[in]    boundary_originx      :    [BOUNDARY] originx
  !! @param[in]    boundary_originy      :    [BOUNDARY] originy
  !! @param[in]    boundary_originz      :    [BOUNDARY] originz
  !! @param[in]    boundary_domain_x     :    [BOUNDARY] domain_x
  !! @param[in]    boundary_domain_y     :    [BOUNDARY] domain_y
  !! @param[in]    boundary_domain_z     :    [BOUNDARY] domain_z
  !! @param[in]    selection_group       :   [SELECTION] group000
  !! @param[in]    restraint_func        :  [RESTRAINTS] func000
  !! @param[in]    restraint_const       :  [RESTRAINTS] const000
  !! @param[in]    restraint_index       :  [RESTRAINTS] index000
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_check_compatible(energy_pairlistdist,   &
                                  energy_table,          &
                                  ensemble_type,         &
                                  boundary_boxsizex,     &
                                  boundary_boxsizey,     &
                                  boundary_boxsizez,     &
                                  boundary_originx,      & 
                                  boundary_originy,      & 
                                  boundary_originz,      & 
                                  boundary_domain_x,     &
                                  boundary_domain_y,     &
                                  boundary_domain_z,     &
                                  selection_group,       &
                                  restraint_func,        &
                                  restraint_const,       &
                                  restraint_index)

    ! formal arguments
    real(wp),                     intent(in)   :: energy_pairlistdist
    logical,                      intent(in)   :: energy_table
    integer,                      intent(in)   :: ensemble_type
    real(wp),                     intent(in)   :: boundary_boxsizex
    real(wp),                     intent(in)   :: boundary_boxsizey
    real(wp),                     intent(in)   :: boundary_boxsizez
    real(wp),                     intent(in)   :: boundary_originx
    real(wp),                     intent(in)   :: boundary_originy
    real(wp),                     intent(in)   :: boundary_originz
    integer,                      intent(in)   :: boundary_domain_x
    integer,                      intent(in)   :: boundary_domain_y
    integer,                      intent(in)   :: boundary_domain_z
    character(*),    allocatable, intent(in)   :: selection_group(:)
    integer,         allocatable, intent(in)   :: restraint_func(:)
    character(*),    allocatable, intent(in)   :: restraint_const(:)
    character(*),    allocatable, intent(in)   :: restraint_index(:)
    
    ! local variables
    integer                :: ip, i, s1, s2
    character(200)         :: file1, file2
    character(50)          :: valstr
    logical                :: lerr, lwrn


    if (.not. main_rank) &
      return

    lerr = .false.
    lwrn = .false.

    if (pio_t0_info%energy_pairlistdist /= energy_pairlistdist) then

      write(MsgOut,'(a,f8.2)') &
           'Pio_Check_Compatible> ERROR : [Energy] pairlist dist: must be ', &
           pio_t0_info%energy_pairlistdist
      lerr = .true.

    end if

    if (pio_t0_info%energy_table .and. .not. energy_table .or. &
        .not. pio_t0_info%energy_table .and. energy_table) then

      if (pio_t0_info%energy_table) then
        valstr = 'YES'
      else
        valstr = 'NO'
      end if

      write(MsgOut,'(a,a)') &
           'Pio_Check_Compatible> ERROR : [Energy] table:         must be ', &
           trim(valstr)
      lerr = .true.

    end if

    if (pio_t0_info%ensemble_type /= ensemble_type) then

      write(MsgOut,'(a,a)') &
           'Pio_Check_Compatible> ERROR : [Ensemble]ensemble: must be ', &
           EnsembleTypes(pio_t0_info%ensemble_type)
      lerr = .true.

    end if

    if (pio_t0_info%boundary_boxsizex /= boundary_boxsizex .or. &
        pio_t0_info%boundary_boxsizey /= boundary_boxsizey .or. &
        pio_t0_info%boundary_boxsizez /= boundary_boxsizez) then

      write(MsgOut,'(a,3f8.2)') &
           'Pio_Check_Compatible> WARNING : [Boundary] box size:    t=0 : ', &
           pio_t0_info%boundary_boxsizex, &
           pio_t0_info%boundary_boxsizey, &
           pio_t0_info%boundary_boxsizez
      lwrn = .true.

    end if

    if (pio_t0_info%boundary_originx /= boundary_originx .or. &
        pio_t0_info%boundary_originy /= boundary_originy .or. &
        pio_t0_info%boundary_originz /= boundary_originz) then

      write(MsgOut,'(a,3f8.2)') &
           'Pio_Check_Compatible> WARNING : [Boundary] origin:      t=0 : ', &
           pio_t0_info%boundary_originx, &
           pio_t0_info%boundary_originy, &
           pio_t0_info%boundary_originz
      lwrn = .true.

    end if

    if (boundary_domain_x == 0 .and. &
        boundary_domain_y == 0 .and. &
        boundary_domain_z == 0) then

      if (pio_t0_info%boundary_domain_xyz /= nproc_world) then
        write(MsgOut,'(a,i8)') &
           'Pio_Check_Compatible> ERROR : # of process :          must be ', &
           pio_t0_info%boundary_domain_xyz
        lerr = .true.
      end if

    else

      if (pio_t0_info%boundary_domain_xyz /= boundary_domain_x * &
                                             boundary_domain_y * &
                                             boundary_domain_z) then
        write(MsgOut,'(a,i8)') &
           'Pio_Check_Compatible> ERROR : [Boundary]domain xyz : must be ',&
           pio_t0_info%boundary_domain_xyz
        lerr = .true.
      end if

    end if

    s1 = 0
    s2 = 0
    if (allocated(pio_t0_info%selection_group)) &
      s1 = size(pio_t0_info%selection_group)
    if (allocated(selection_group)) &
      s2 = size(selection_group)

    if (s1 /= s2) then
      
      write(MsgOut,'(a,i8)') &
           'Pio_Check_Compatible> ERROR : [Selection] # of groups : must be ',&
           s1
      lerr = .true.

    else

      do i = 1, s1
        if (pio_t0_info%selection_group(i) /= selection_group(i)) then
          write(MsgOut,'(a,i8,a,a)') &
               'Pio_Check_Compatible> ERROR : [Selection] group ', i, &
               ' : must be ', trim(pio_t0_info%selection_group(i))
          lerr = .true.

        end if
      end do

    end if

    s1 = 0
    s2 = 0
    if (allocated(pio_t0_info%restraint_func)) &
      s1 = size(pio_t0_info%restraint_func)
    if (allocated(restraint_func)) &
      s2 = size(restraint_func)

    if (s1 /= s2) then
      
      write(MsgOut,'(a,i8)') &
           'Pio_Check_Compatible> ERROR : [Restraint] # of funcs : must be ',&
           s1
      lerr = .true.

    else

      do i = 1, s1
        if (pio_t0_info%restraint_func(i) /= restraint_func(i)) then
          write(MsgOut,'(a,i8,a,i8)') &
               'Pio_Check_Compatible> ERROR : [Restraint] func ', i, &
               ' : must be ', pio_t0_info%restraint_func(i)
          lerr = .true.

        end if
      end do

      do i = 1, s1
        if (pio_t0_info%restraint_const(i) /= restraint_const(i)) then
          write(MsgOut,'(a,i8,a,a)') &
               'Pio_Check_Compatible> ERROR : [Restraint] const ', i, &
               ' : must be ', trim(pio_t0_info%restraint_const(i))
          lerr = .true.

        end if
      end do

      do i = 1, s1
        if (pio_t0_info%restraint_index(i) /= restraint_index(i)) then
          write(MsgOut,'(a,i8,a,a)') &
               'Pio_Check_Compatible> ERROR : [Restraint] index ', i, &
               ' : must be ', trim(pio_t0_info%restraint_index(i))
          lerr = .true.

        end if
      end do

    end if


    if (lerr) &
      call error_msg('Pio_Check_Compatible> ERROR')

    if (lwrn) &
      write(MsgOut,*) ' '

    return

  end subroutine pio_check_compatible

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_check_ranked_file
  !> @brief        check ranked file name or not
  !! @authors      NT
  !! @param[in]    filename : file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_check_ranked_file(filename)

    ! function
    logical                  :: pio_check_ranked_file

    ! formal arguments
    character(100),          intent(in) :: filename

    ! local variables
    integer                  :: ci1, ci2


    ! check filename
    !
    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    pio_check_ranked_file = (ci1 /= 0 .and. ci2 /= 0 .and. ci1 < ci2)

    return

  end function pio_check_ranked_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_get_ranked_filename
  !> @brief        get ranked file name
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[out]   nplace   : number of places
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_get_ranked_filename(filename, nplace)

    ! function
    character(100)           :: pio_get_ranked_filename

    ! formal arguments
    character(*),            intent(in) :: filename
    integer,       optional, intent(in) :: nplace

    ! local variables
    integer                  :: ci1, ci2, ip
    character(50)            :: fmt_str

    ! constants
    integer,       parameter :: Places (7) = &
         (/10, 100, 1000, 10000, 100000, 1000000, 10000000/)


    if (present(nplace)) then
      ip = nplace
    else
      do ip = 1, size(Places)
        if (nproc_world < Places(ip)) &
          exit
      end do
    end if

    ! check filename
    !
    pio_get_ranked_filename = filename
    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    if (ci1 == 0 .or. ci2 ==0 .or. ci1 > ci2) then
      call error_msg( &
      'Pio_Get_Ranked_Filename> Filename is not correctly ranked in [OUTPUT]')
    end if

    write(fmt_str,'(A,I0,A,I0,A)') '(A,I',ip,'.',ip,',A)'

    write(pio_get_ranked_filename,fmt=fmt_str)   &
         trim(filename(1:ci1-1)), my_world_rank, &
         trim(filename(ci2+1:))

    return

  end function pio_get_ranked_filename

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pio_open_file
  !> @brief        open domain decomposition restart file
  !! @authors      NT
  !! @param[inout] file     : file unit number
  !! @param[in]    filename : file name
  !! @param[in]    in_out   : input or output
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pio_open_file(file, filename, in_out)

    ! formal arguments
    integer,                 intent(inout) :: file
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: in_out

    ! local variables
    integer                  :: imark, ist


    select case(in_out)

    case (IOFileInput)

      call open_binary_file(file, filename, in_out, &
                            pio_get_file_endian(filename))
      read(file, iostat=ist) imark

      if (ist /= 0 .or. imark /= PioFileHeaderMarker) &
        call error_msg('Pio_Open_File> Unknown file format :'//&
        trim(filename))

    case (IOFileOutputNew, IOFileOutputReplace)

      call open_binary_file(file, filename, in_out)
      write(file) PioFileHeaderMarker

    case default

      call error_msg('Pio_Open_File> Unknown I/O mode')

    end select

    return

  end subroutine pio_open_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_get_file_endian
  !> @brief      
  !! @authors      NT
  !! @param[in]  
  !! @param[out] 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_get_file_endian(filename)

    ! function
    integer          :: pio_get_file_endian

    ! formal arguments
    character(*),    intent(in) :: filename

    ! local variables
    integer          :: file, ival


    file = get_unit_no()

    open(file, &
         file    = filename,        &
         status  = 'old',           &
         form    = 'unformatted',   &
         convert = 'little_endian', &
         access  = 'direct',        &
         recl    = 4)
    read(file,rec=1) ival

    if (ival /= 4) then

      close(file)

      open(file, &
           file    = filename,      &
           status  = 'old',         &
           form    = 'unformatted', &
           convert = 'big_endian',  &
           access  = 'direct',      &
           recl    = 4)
      read(file,rec=1) ival

      if (ival /= 4) &
        call error_msg('Pio_Get_File_Endian> unknown file format.')

      pio_get_file_endian = IOFileBigEndian

    else

      pio_get_file_endian = IOFileLittleEndian

    end if

    close(file)
    call free_unit_no(file)

    return

  end function pio_get_file_endian

end module cg_parallel_io_mod
