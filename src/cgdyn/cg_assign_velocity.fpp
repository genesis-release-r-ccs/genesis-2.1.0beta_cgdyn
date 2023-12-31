!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_assign_velocity_mod
!> @brief   generate initial velocities and remove trans- and rotational motions
!! @authors Jaewoon Jung (JJ)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_assign_velocity_mod

  use random_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use cg_dynvars_mod
  use cg_domain_str_mod
  use cg_dynamics_str_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! subroutines
  public  :: initial_velocity
  public  :: stop_trans_rotation
  private :: reduce_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    initial_velocity
  !> @brief        generate initial velocities
  !! @authors      JJ
  !! @param[in]    firstt   : initial temperature
  !! @param[in]    natom    : number of atom
  !! @param[in]    ncell    : total number of cells in each domain
  !! @param[in]    id_g2l   : map of atom index (global -> local)
  !! @param[in]    mass     : mass in each domain
  !! @param[inout] iseed    : random number seed
  !! @param[inout] velocity : velocities generated
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine initial_velocity(firstt, iseed, domain)

    ! formal arguments
    real(wip),               intent(in)    :: firstt
    integer,                 intent(inout) :: iseed
    type(s_domain), target,  intent(inout) :: domain

    ! pointer
    integer,         pointer :: id_g2l(:)
    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: velocity(:,:)

    ! local variables
    real(wip)                :: alpha1, beta1, alpha2, beta2, alpha3, beta3
    real(wip)                :: sigma, myu, kBT
    integer                  :: i, j
    integer                  :: natom, natom_domain

    id_g2l       => domain%id_g2l
    mass         => domain%mass
    velocity     => domain%velocity
    natom        =  domain%num_atom_all
    natom_domain =  domain%num_atom_domain

    if (main_rank .or. replica_main_rank) then
      write(DynvarsOut,'(A)') 'Initial_Velocity> Generate initial velocities'
      write(DynvarsOut,'(A20,I10)')   '  iseed           = ', iseed
      write(DynvarsOut,'(A20,F10.3)') '  temperature     = ', firstt
      write(DynvarsOut,'(A)') ' '
    end if

    ! generate initial velocities 
    !
    kBT = KBOLTZ * firstt
    myu = 0.0_wip

    do i = 1, natom

       alpha1 = sqrt( -2.0_wip * log(random_get_legacy(iseed)) )
       beta1  = cos ( 2.0_wip * PI * random_get_legacy(iseed) )

       alpha2 = sqrt( -2.0_wip * log(random_get_legacy(iseed)) )
       beta2  = cos ( 2.0_wip * PI * random_get_legacy(iseed) )

       alpha3 = sqrt( -2.0_wip * log(random_get_legacy(iseed)) )
       beta3  = cos ( 2.0_wip * PI * random_get_legacy(iseed) )

       j = id_g2l(i)
       if (j <= natom_domain .and. j > 0) then
          if (abs(mass(j)) > EPS) then
            sigma = sqrt( kBT / mass(j) )
          else
            sigma = 0.0_wip
          end if
          velocity(j,1) = sigma * alpha1 * beta1 + myu
          velocity(j,2) = sigma * alpha2 * beta2 + myu
          velocity(j,3) = sigma * alpha3 * beta3 + myu

       end if

    end do
  
    return 

  end subroutine initial_velocity

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    stop_trans_rotation
  !> @brief        remove trans and rotational motion about the center of mass
  !! @authors      JJ
  !! @param[in]    ncell      : number of cell
  !! @param[in]    natom      : number of atom in each cell
  !! @param[in]    stop_trans : flag for stop translational motion
  !! @param[in]    stop_rot   : flag for stop rotational motion
  !! @param[in]    mass       : mass in each domain
  !! @param[in]    coord      : coordinates in each domain
  !! @param[inout] velocity   : velocities in each domain
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine stop_trans_rotation(dynamics, domain)

    ! formal arguments
    type(s_dynamics),         intent(inout) :: dynamics
    type(s_domain),   target, intent(inout) :: domain

    ! pointer
    real(wip),       pointer :: mass(:)
    real(wip),       pointer :: coord(:,:), velocity(:,:)

    ! local variables
    real(dp)                 :: inertia(3,3), inv_inertia(3,3)
    real(dp)                 :: total_mass, ccom(3), vcom(3)
    real(dp)                 :: l_angular(1:3)
    real(dp)                 :: c(3), v(3), omega(3)
    real(dp)                 :: xx, xy, xz, yy, yz, zz
    integer                  :: i
    integer                  :: num_atom_domain
    logical                  :: stop_trans, stop_rot

    coord           => domain%coord
    mass            => domain%mass
    velocity        => domain%velocity
    num_atom_domain =  domain%num_atom_domain
    stop_trans      =  dynamics%stop_com_translation
    stop_rot        =  dynamics%stop_com_rotation

    ! calculate the center of mass of coordinates and velocities
    !
    ccom(1:3)  = 0.0_dp
    vcom(1:3)  = 0.0_dp
    total_mass = 0.0_dp
    do i = 1, num_atom_domain
      ccom(1)  = ccom(1) + coord(i,1)*mass(i)
      ccom(2)  = ccom(2) + coord(i,2)*mass(i)
      ccom(3)  = ccom(3) + coord(i,3)*mass(i)
      vcom(1)  = vcom(1) + velocity(i,1)*mass(i)
      vcom(2)  = vcom(2) + velocity(i,2)*mass(i)
      vcom(3)  = vcom(3) + velocity(i,3)*mass(i)
      total_mass = total_mass + mass(i)
    end do

    call reduce_com(ccom, vcom, total_mass)
    
    ccom(1:3) = ccom(1:3) / total_mass
    vcom(1:3) = vcom(1:3) / total_mass

    ! calculate the angular momentum
    !
    l_angular(1:3) = 0.0_dp
    do i = 1, num_atom_domain
      c(1) = coord(i,1) - ccom(1)
      c(2) = coord(i,2) - ccom(2)
      c(3) = coord(i,3) - ccom(3)
      v(1) = velocity(i,1) - vcom(1)
      v(2) = velocity(i,2) - vcom(2)
      v(3) = velocity(i,3) - vcom(3)
      l_angular(1) = l_angular(1) + (c(2)*v(3) - c(3)*v(2))*mass(i)
      l_angular(2) = l_angular(2) + (c(3)*v(1) - c(1)*v(3))*mass(i)
      l_angular(3) = l_angular(3) + (c(1)*v(2) - c(2)*v(1))*mass(i)
    end do
 
#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, l_angular, 3, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    ! calculate the inertia tensor (I)
    !
    xx = 0.0_dp
    xy = 0.0_dp
    xz = 0.0_dp
    yy = 0.0_dp
    yz = 0.0_dp
    zz = 0.0_dp

    do i = 1, num_atom_domain
      c(1) = coord(i,1) - ccom(1)
      c(2) = coord(i,2) - ccom(2)
      c(3) = coord(i,3) - ccom(3)
      xx = xx + c(1)*c(1)*mass(i)
      xy = xy + c(1)*c(2)*mass(i)
      xz = xz + c(1)*c(3)*mass(i)
      yy = yy + c(2)*c(2)*mass(i)
      yz = yz + c(2)*c(3)*mass(i)
      zz = zz + c(3)*c(3)*mass(i)
    end do

    inertia(1,1) =   yy + zz
    inertia(1,2) = - xy
    inertia(1,3) = - xz
    inertia(2,1) = - xy
    inertia(2,2) =   xx + zz
    inertia(2,3) = - yz
    inertia(3,1) = - xz
    inertia(3,2) = - yz
    inertia(3,3) =   xx + yy

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(mpi_in_place, inertia, 9, mpi_real8, &
                       mpi_sum, mpi_comm_country, ierror)
#endif

    ! calculate the inverse matrix of the inertia tensor
    !
    call compute_inverse_matrix(3, inertia, inv_inertia)
 
    ! calculate the angular velocity (OMG)
    !
    omega(1) = inv_inertia(1,1)*l_angular(1) + inv_inertia(1,2)*l_angular(2) &
             + inv_inertia(1,3)*l_angular(3)
    omega(2) = inv_inertia(2,1)*l_angular(1) + inv_inertia(2,2)*l_angular(2) &
             + inv_inertia(2,3)*l_angular(3)
    omega(3) = inv_inertia(3,1)*l_angular(1) + inv_inertia(3,2)*l_angular(2) &
             + inv_inertia(3,3)*l_angular(3)

    ! remove translational motion
    !
    if (stop_trans) then
      do i = 1, num_atom_domain
        velocity(i,1) = velocity(i,1) - vcom(1)
        velocity(i,2) = velocity(i,2) - vcom(2)
        velocity(i,3) = velocity(i,3) - vcom(3)
      end do
    end if

    ! remove rotational motion
    !
    if (stop_rot) then
      do i = 1, num_atom_domain
        c(1) = coord(i,1) - ccom(1)
        c(2) = coord(i,2) - ccom(2)
        c(3) = coord(i,3) - ccom(3)
        velocity(i,1) = velocity(i,1) - (omega(2)*c(3) - omega(3)*c(2))
        velocity(i,2) = velocity(i,2) - (omega(3)*c(1) - omega(1)*c(3))
        velocity(i,3) = velocity(i,3) - (omega(1)*c(2) - omega(2)*c(1))
      end do
    end if

    return

  end subroutine stop_trans_rotation

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reduce_com
  !> @brief
  !! @authors      JJ
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reduce_com(val1, val2, val3)

    ! formal arguments
    real(dp),                intent(inout) :: val1(:), val2(:), val3

    ! local variables
    real(dp)                      :: before_reduce(7), after_reduce(7)


    before_reduce(1:3) = val1(1:3)
    before_reduce(4:6) = val2(1:3)
    before_reduce(7)   = val3

#ifdef HAVE_MPI_GENESIS
    call mpi_allreduce(before_reduce, after_reduce, 7, mpi_real8, &
                       mpi_sum, mpi_comm_city, ierror)
#else
    after_reduce(1:7) = before_reduce(1:7)
#endif

    val1(1:3)    = after_reduce(1:3)
    val2(1:3)    = after_reduce(4:6)
    val3         = after_reduce(7)

    return

  end subroutine reduce_com


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    compute_inverse_matrix
  !> @brief        Gauss-Jordan elimination method with partial pivoting
  !! @authors      TM
  !! @param[in]    n     : matrix dimension
  !! @param[in]    m     : original matrix
  !! @param[out]   inv_m : inverse matrix
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine compute_inverse_matrix(n, m, inv_m)

    ! formal arguments
    integer,      intent(in)  :: n
    real(dp),     intent(in)  :: m(n,n)
    real(dp),     intent(out) :: inv_m(n,n)

    ! local variables
    integer                   :: i, j, k, imax
    real(dp)                  :: a(n,2*n), amax, swap, coef, pivot

    do i = 1, n
      do j = 1, n
        a(i,j) = m(i,j)
        if (i == j) a(i,j+n) = 1.0_dp
        if (i /= j) a(i,j+n) = 0.0_dp
      end do
    end do

    do k = 1, n
      amax = a(k,k)
      imax = k
      do i = 1, n
        if (abs(a(i,k)) > amax) then
          amax = a(i,k)
          imax = i
        end if
      end do

      do j = k, 2*n
        swap = a(k,j)
        a(k,j) = a(imax,j)
        a(imax,j) = swap
      end do

      pivot = 1.0_dp/a(k,k)
      do j = k, 2*n
        a(k,j) = a(k,j)*pivot
      end do
      do i = 1, n
        if (i /= k) then
          coef = a(i,k)
          do j = k, 2*n
            a(i,j) = a(i,j) - a(k,j)*coef
          end do
        end if
      end do
    end do

    do i = 1, n
      do j = 1, n
        inv_m(i,j) = a(i,j+n)
      end do
    end do

    return

  end subroutine compute_inverse_matrix

end module cg_assign_velocity_mod
