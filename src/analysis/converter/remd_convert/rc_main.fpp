!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  rc_main
!! @brief   convert REMD trajectory files
!! @authors Norio Takase (NT), Yuji Sugita (YS)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

program rc_main

  use rc_convert_mod
  use rc_setup_mod
  use rc_control_mod
  use rc_option_str_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_str_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none

  ! local variables
  character(MaxFilename) :: ctrl_filename
  type(s_ctrl_data)      :: ctrl_data
  type(s_molecule)       :: molecule
  type(s_trajectory)     :: trajectory
  type(s_fitting)        :: fitting
  type(s_option)         :: option
  type(s_input)          :: input
  type(s_output)         :: output


  my_city_rank = 0
  nproc_city   = 1
  main_rank    = .true.


  ! show usage
  !
  call usage(ctrl_filename)


  ! [Step1] Read control file
  !
  write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Convert'
  write(MsgOut,'(A)') ' '

  call control(ctrl_filename, ctrl_data)


  ! [Step2] Set relevant variables and structures 
  !
  write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
  write(MsgOut,'(A)') ' '

  call setup(ctrl_data,  &
             molecule,   &
             trajectory, &
             fitting,    &
             option,     &
             input,      &
             output)


  ! [Step3] Convert REMD trajectory files
  !
  write(MsgOut,'(A)') '[STEP3] Convert REMD trajectory files'
  write(MsgOut,'(A)') ' '

  call convert(input,      &
               output,     &
               option,     &
               molecule,   &
               trajectory, &
               fitting)


  ! [Step4] Deallocate memory
  !
  write(MsgOut,'(A)') '[STEP4] Deallocate memory'
  write(MsgOut,'(A)') ' '

  call dealloc_option(option)
  call dealloc_fitting(fitting)
  call dealloc_trajectory(trajectory)
  call dealloc_molecules_all(molecule)

  stop

end program rc_main