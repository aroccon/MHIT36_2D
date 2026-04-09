program interface_indicator
  implicit none

  integer,  parameter :: nx = 2048, ny = 2048
  real(8),  parameter :: lx = 2048.d0, ly = 2048.d0
  real(8),  parameter :: dx = lx / dble(nx)
  real(8),  parameter :: dy = ly / dble(ny)
  integer,  parameter :: t_start = 0
  integer,  parameter :: t_end   = 10000000
  integer,  parameter :: t_step  = 4000
  real(8), allocatable :: phi_field(:,:)
  integer :: i, j, snap, im, ip, jm, jp
  real(8) :: dphidx, dphidy, I_sum, I_val
  character(len=64) :: fname_phi
  allocate( phi_field(nx, 0:ny+1) )

  open(unit=30, file='interface_timeseries.dat', status='replace')
  write(30,'(A)') '# snap   I'

  do snap = t_start, t_end, t_step

    write(fname_phi, '("phi_", I8.8, ".dat")') snap

    open(unit=10, file='../../src/output/' // trim(fname_phi), form='unformatted', access='stream', status='old')
    read(10) phi_field; close(10)


    I_sum = 0.d0
    do j = 1, ny
      do i = 1, nx
        im = i - 1
        ip = i + 1
        jm = j - 1
        jp = j + 1
        if (im < 1) im = 1
        if (ip > nx) ip = nx
        dphidx = (phi_field(ip,j) - phi_field(im,j)) / (2.d0*dx)
        dphidy = (phi_field(i,jp) - phi_field(i,jm)) / (2.d0*dy)
        I_sum  = I_sum + dsqrt(dphidx**2 + dphidy**2)
      end do
    end do

    I_val = I_sum / dble(nx*ny)

    write(30, '(I10, 2X, ES15.7)') snap, I_val
    write(*,  '("Snapshot ",I8," I=",ES15.7)') snap, I_val

  end do

  close(30)
  deallocate(phi_field)

end program interface_indicator