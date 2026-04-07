program nusselt_rb
  implicit none

  integer,  parameter :: nx = 2048, ny = 2048
  real(8),  parameter :: ly = 2048.d0
  real(8),  parameter :: kappa = 0.166d0
  real(8),  parameter :: deltaT = 1.d0   ! T_hot - T_cold
  integer,  parameter :: t_start = 0
  integer,  parameter :: t_end   = 10000000
  integer,  parameter :: t_step  = 4000
  real(8), allocatable :: t_field(:,:), v_field(:,:)
  integer :: i, j, snap, n_snaps
  real(8) :: v_center, vT_sum, nu, nu_run_avg
  character(len=64) :: fname_t, fname_v
  allocate( t_field(nx, 0:ny+1) )
  allocate( v_field(nx,   ny+1) )

  n_snaps    = 0

  open(unit=30, file='nusselt_timeseries.dat', status='replace')
  write(30,'(A)') '# snap   Nu'

  do snap = t_start, t_end, t_step

    ! read fields
    write(fname_t, '("t_", I8.8, ".dat")') snap
    write(fname_v, '("v_", I8.8, ".dat")') snap

    open(unit=10, file='../../src/output/' // trim(fname_t), form='unformatted', access='stream', status='old')
    read(10) t_field; close(10)
    open(unit=11, file='../../src/output/' // trim(fname_v), form='unformatted', access='stream', status='old')
    read(11) v_field; close(11)

    ! compute <vT> over full domain
    vT_sum = 0.d0
    do j = 1, ny
      do i = 1, nx
        v_center = 0.5d0 * (v_field(i,j) + v_field(i,j+1))
        vT_sum   = vT_sum + v_center * t_field(i,j)
      end do
    end do

    nu = 1.d0 + (vT_sum / dble(nx*ny))/(kappa*deltaT/ly)

    ! write to time series
    write(30, '(I10, 1(2X, ES15.7))') snap, nu
    write(*,'("Processed snapshot ",I8)') snap
  end do
  close(30)

  deallocate(t_field, v_field)

end program nusselt_rb