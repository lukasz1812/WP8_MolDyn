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
  write(*,*) class, "",atom, mass, natom, length, steps, stept



!Case selection for different systems

select case (class)

  case ("SC")
    call SC_Grid(natom,length,atom)

  case ("sc")
    call SC_Grid(natom,length,atom)

  case ("PC")
    call SC_Grid(natom,length,atom)

  case ("pc")
    call SC_Grid(natom,length,atom)

  case ("fcc")
    call FCC_Grid(natom,length,atom)

  case ("FCC")
    call FCC_Grid(natom,length,atom)

end select



end program MolDyn

!-------------------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@---Subroutines---@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------------------



!=====================================Subroutine for reading the input datas from the file=====================================
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


!=====================================Subroutine for the primitive cubic =====================================
subroutine SC_Grid(natom,length,atom)

  !Declaration of local Variables
  integer :: natom,n,nl,i,j,k
  real ::  hl,dl,x,y,z, length
  real, allocatable :: rx(:),ry(:),rz(:)
  character :: atom*2
  allocate(rx(natom),ry(natom),rz(natom))
  write(*,*) "length in subroutine", length

  hl=length/2
  nl=int(natom**(1.0d0/3.0d0))


if (nl**3.lt.natom) then
    nl=nl+1
end if
write (*,*) nl

dl=length/nl

n=0
  do i=0, nl-1

    do j=0, nl-1

      do k=0, nl-1

          n=n+1
          if(natom>=n) then

            x=i*dl-hl
            rx(n)=x
            y=j*dl-hl
            ry(n)=y
            z=k*dl-hl
            rz(n)=z
            write(14,*) atom, x, y, z

          end if

        end do

      end do

    end do

end subroutine SC_Grid
!---------------------------- --End of primitive cubic grid subroutine------------------------------


!=====================================Subroutine for the face centered cubic grid=====================================
subroutine FCC_Grid(natom,length,atom)

  !Declaration of local Variables
  integer :: natom,n,nl,i,j,k,sum
  real ::  hl,dl,x,y,z, length
  real, allocatable :: rx(:),ry(:),rz(:)
  character :: atom*2
  allocate(rx(natom),ry(natom),rz(natom))

  write(*,*) "length in subroutine", length


hl=length/2

nl=int((0.25*natom)**(1.0d0/3.0d0))

write (*,*) nl, "nl1"


  if(4*(nl**3)<natom) then

    nl =nl+1

  end if
  write (*,*) nl, "nl2"
  nl=nl


dl=length/nl
write(*,*) dl

n=0




do i=0,2*nl-2
  do j=0,2*nl-2
    do k=0,2*nl-2


      sum=i+j+k

      if (modulo(sum,2)==0) then !checking if the sum of all coordinates is even

        n=n+1
          if (natom>=n) then




            x=i*dl-hl
            rx(n)=x
            y=j*dl-hl
            ry(n)=y
            z=k*dl-hl
            rz(n)=z
            write(14,*) atom,x,y,z, i, j, k   !prints out “fort.14” file with coordinates

            end if
          end if
        end do
      end do
    end do



end subroutine FCC_grid
!---------------------------- --End of face ceneterd cubic grid subroutine------------------------------


!=====================================Force calculation subroutine=====================================
subroutine calc_force(natom,rx,ry,rz,fx,fy,fz,boxl,ax,ay,az)

  !Declaration of local Variables
  integer :: i, natom
  real, allocatable :: fx(:), fy(:), fz(:)
  real, allocatable ::rx(:), ry(:), rz(:), ax(:), ay(:), az(:)
  real :: k


  allocate(fx(natom), fy(natom), fz(natom))
  allocate(rx(natom), ry(natom), rz(natom), ax(natom), ay(natom), az(natom))

    do i=1,natom !Calculating forces for each of natom and writing it into correspondig vector element

        fx(i)=0.0d0;fy(i)=0.0d0;fz(i)=0.0d0
        fx(i)=-k*(rx(i)-ax(i))
        fy(i)=-k*(ry(i)-ay(i))
        fz(i)=-k*(rz(i)-az(i))

    end do

    !deallocation of memory
    deallocate(fx,fy,fz)
    deallocate(rx,ry,rz,ax,ay,az)

  end subroutine calc_force
  !---------------------------- --End of force calculation subroutine------------------------------


!=====================================POtential energy calculation subroutine=====================================
subroutine calc_harmpot(natom,rx,ry,rz,v,boxl,ax,ay,az,pot_harm)

    !Declaration of local Variables
    integer :: i, natom
    real, allocatable ::rx(:), ry(:), rz(:), ax(:), ay(:), az(:)
    real :: pot_harm

    allocate(rx(natom), ry(natom), rz(natom), ax(natom), ay(natom), az(natom))

    pot_harm=0


      do i=1,natom !Calculating potential energy for every atom

          pot_harm=pot_harm+0.50d0*((rx(i)-ax(i)**2+ry(i)-ay(i)**2+rz(i)-az(i)**2))

      end do

      !deallocation of memory
      deallocate(rx,ry,rz,ax,ay,az)

    end subroutine calc_harmpot
!------------------------------End of potential energy calculation subroutine------------------------------


!=====================================Velocity Verlet Subroutine=====================================
subroutine vverlet(natom)

  !declaration of local Variables
  integer :: i, natom,steps, k
  real :: kin,mass, delta, etot, pot_harm
  real, allocatable ::rx(:), ry(:), rz(:), vx(:), vy(:), vz(:)

  allocate(rx(natom), ry(natom), rz(natom), vx(natom), vy(natom), vz(natom))

  call calc_force(natom,rx,ry,rz,fx,fy,fz,boxl,ax,ay,az)
  call calc_harmpot(natom,rx,ry,rz,v,boxl,ax,ay,az, pot_harm)

    do i=1, natom !Setting velocities for every atom in every dimension as 0

        vx(i)=0.0d0;vy(i)=0.0d0;vz(i)=0.0d0

    end do

    kin=0.0d0
    k=0.0d0

    do i=1,natom !Summing the kinetical energy over all atoms

      kin=kin+0.5*mass*(vx(i)**2+vy(i)**2+vz(i)**2)

    end do
    
open(unit=file_no1,file="kinetic.txt",action="write")
open(unit=file_no2,file="potential_harmonic.txt",action="write")
open(unit=file_no3,file="energy_total.txt",action="write")
  
    do k=1, steps !Main loop

        kin=0.0d0

          do i=1,natom

            !Calculating new velocities for every atom in every dimension
            vx(i)=vx(i)+0.5*(fx(i)/mass)*delta
            vy(i)=vy(i)+0.5*(fy(i)/mass)*delta
            vz(i)=vz(i)+0.5*(fz(i)/mass)*delta

            !Calculating new positions for every atom in every dimension
            rx(i)=rx(i)+vx(i)*delta
            ry(i)=ry(i)+vy(i)*delta
            rz(i)=rz(i)+vz(i)*delta

          end do

          call calc_force(natom,rx,ry,rz,fx,fy,fz,boxl,ax,ay,az)
            

          do i=1, natom !rescale velocities

            vx(i)=vx(i)+0.5*(fx(i)/mass)*delta
            vy(i)=vy(i)+0.5*(fy(i)/mass)*delta
            vz(i)=vz(i)+0.5*(fz(i)/mass)*delta

            !calculate sum over all kinetical energies
            kin=kin+0.5*mass*(vx(i)**2+vy(i)**2+vz(i)**2)
            write(file_no1,*) "Iteration", k, "Energy", kin

          end do

            call calc_harmpot(natom,rx,ry,rz,v,boxl,ax,ay,az,pot_harm)
            write(file_no2,*) "Iteration", k, "Energy", pot_harm

            etot=pot_harm+kin
            write(file_no3,*) "Iteration", k, "Energy", etot
            
            
    end do
    close(file_no1)
    close(file_no2)
    close(file_no3)

end subroutine vverlet
!------------------------------End of velocity Verlet subroutine------------------------------
