!This programm  was written by Łukasz Wantoch
!as a part of internship of module WP8
!Molecular Dynamics of Time Dependent Phenomena
!@ the Bonn University in Winterterm 2020/21


program MolDyn
  use iso_fortran_env, only : output_unit, error_unit

  !Declaration of global variables
  character :: filename*100


  !Small user interface
  write(*,*) "Write te input data File name:"
  read(*,*) filename

end program MolDyn



!-------------------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---Subroutines---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------------------



!Subroutine for reading the input datas from the file
subroutine reader(filename,class,atom,mass,natom,length,steps,stept)

  !Declaration of local variables


end subroutine reader
