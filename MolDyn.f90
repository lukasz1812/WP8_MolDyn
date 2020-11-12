!This programm  was written by Łukasz Wantoch
!as a part of internship of module WP8
!Molecular Dynamics of Time Dependent Phenomena
!@ the Bonn University in Winterterm 2020/21


program MolDyn
  use iso_fortran_env, only : output_unit, error_unit

  !Declaration of global variables
  character :: filename*100,class*3,atom*2
  integer :: natom, steps
  real :: mass, length, stept

  !Small user interface
  write(*,*) "Write te input data File name:"
  read(*,*) filename
  call input_reader(filename,class,atom,mass,natom,length,steps,stept)

  !Start variable check
  write(*,*) class, atom, mass, natom, length, steps, stept
  
end program MolDyn



!-------------------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---Subroutines---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------------------



!Subroutine for reading the input datas from the file
subroutine input_reader(filename,class,atom,mass,natom,length,steps,stept)

  !Declaration of local variables
  character :: filename*100,class*3,atom*2
  integer :: natom, steps, io
  real :: mass, length, stept

  !open input file
  open(file=fileName,newunit=io)
  
  !Input of start Variables
  read(io,*) class, atom, mass, natom, length, steps, stept
  close(io)

end subroutine input_reader
!---------------------------- --End of input reader subroutine------------------------------
