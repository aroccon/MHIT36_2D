program para
! Convert phi_********.dat and t_********.dat snapshots (unformatted stream,
! written by src/readwrite.f90) into legacy VTK STRUCTURED_POINTS files that
! ParaView can load directly, one per snapshot, with phi and temperature as
! point scalars on the cartesian cell-centered grid.
!
! Grid/domain parameters below must match src/module.f90 (nx,ny) and
! src/input.inp (lx,ly) for the run being post-processed.
implicit none

integer,  parameter :: nx = 2048, ny = 512
real(8),  parameter :: lx = 4096.d0, ly = 1024.d0
real(8),  parameter :: dx = lx/dble(nx), dy = ly/dble(ny)

! snapshot range/stride to scan (must match dump used in src/input.inp)
integer,  parameter :: t_start = 0
integer,  parameter :: t_end   = 4000000
integer,  parameter :: t_step  = 8000

character(len=*), parameter :: datdir = '../../src/output/'
character(len=*), parameter :: outdir = 'vtk/'

real(8), allocatable :: phi_field(:,:), t_field(:,:)
integer :: snap
character(len=128) :: fname_phi, fname_t, fname_out
logical :: has_phi, has_t

allocate(phi_field(nx,0:ny+1))
allocate(t_field(nx,0:ny+1))

call execute_command_line('mkdir -p '//outdir)

do snap = t_start, t_end, t_step

  write(fname_phi,'(a,"phi_",i8.8,".dat")') datdir, snap
  write(fname_t,  '(a,"t_",  i8.8,".dat")') datdir, snap

  inquire(file=trim(fname_phi), exist=has_phi)
  inquire(file=trim(fname_t),   exist=has_t)
  if (.not. (has_phi .and. has_t)) cycle

  open(unit=10, file=trim(fname_phi), form='unformatted', access='stream', status='old', convert='little_endian')
  read(10) phi_field
  close(10)

  open(unit=11, file=trim(fname_t), form='unformatted', access='stream', status='old', convert='little_endian')
  read(11) t_field
  close(11)

  write(fname_out,'(a,"field_",i8.8,".vtk")') outdir, snap
  call write_vtk(trim(fname_out), nx, ny, dx, dy, phi_field, t_field)

  write(*,'("Written snapshot ",i8," -> ",a)') snap, trim(fname_out)

end do

deallocate(phi_field, t_field)

contains

subroutine write_vtk(fname, nx, ny, dx, dy, phi, temp)
! Legacy VTK binary STRUCTURED_POINTS file (big-endian, as required by the format)
implicit none
character(len=*), intent(in) :: fname
integer,  intent(in) :: nx, ny
real(8),  intent(in) :: dx, dy
real(8),  intent(in) :: phi(nx,0:ny+1), temp(nx,0:ny+1)
integer :: i, j, u
character(len=256) :: line
real(4) :: val

u = 20
open(unit=u, file=fname, form='unformatted', access='stream', status='replace', convert='big_endian')

write(u) '# vtk DataFile Version 3.0'//char(10)
write(u) 'MHIT36_2D phase field and temperature'//char(10)
write(u) 'BINARY'//char(10)
write(u) 'DATASET STRUCTURED_POINTS'//char(10)

write(line,'("DIMENSIONS ",i0," ",i0," 1")') nx, ny
write(u) trim(line)//char(10)

write(line,'("ORIGIN ",es15.7," ",es15.7," 0.0")') dx/2.d0, dy/2.d0
write(u) trim(line)//char(10)

write(line,'("SPACING ",es15.7," ",es15.7," 1.0")') dx, dy
write(u) trim(line)//char(10)

write(line,'("POINT_DATA ",i0)') nx*ny
write(u) trim(line)//char(10)

write(u) 'SCALARS phi float 1'//char(10)
write(u) 'LOOKUP_TABLE default'//char(10)
do j = 1, ny
  do i = 1, nx
    val = real(phi(i,j),4)
    write(u) val
  end do
end do
write(u) char(10)

write(u) 'SCALARS temperature float 1'//char(10)
write(u) 'LOOKUP_TABLE default'//char(10)
do j = 1, ny
  do i = 1, nx
    val = real(temp(i,j),4)
    write(u) val
  end do
end do
write(u) char(10)

close(u)

end subroutine write_vtk

end program para
