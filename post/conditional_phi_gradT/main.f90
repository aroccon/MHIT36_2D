program cond_avg_phi_gradT
  implicit none
  integer,  parameter :: nx = 2048, ny = 2048
  integer,  parameter :: t_start = 2000000
  integer,  parameter :: t_end   = 10000000
  integer,  parameter :: t_step  = 4000
  integer,  parameter :: n_bins  = 100

  real(8),  parameter :: Lx = 1.d0   ! domain length in x -- adjust
  real(8),  parameter :: Ly = 1.d0   ! domain length in y -- adjust
  real(8),  parameter :: dx = Lx / dble(nx)
  real(8),  parameter :: dy = Ly / dble(ny)

  real(8), allocatable :: t_field(:,:), phi_field(:,:)
  real(8), allocatable :: cond_phi_sum(:), cond_count(:)

  integer :: i, j, snap, ibin, n_snaps
  real(8) :: dTdx, dTdy, gradT_mag
  real(8) :: bin_width, gradT_bin_center, phi_cond
  real(8) :: gradT_max
  character(len=128) :: fname_t, fname_phi

  ! ----------------------------------------------------------------
  ! gradT_max: set this to the maximum |gradT| you expect.
  ! A good estimate is (t_hot - t_cold) / thickness_of_BL
  ! or simply run once and check the output range.
  ! ----------------------------------------------------------------
  gradT_max = 100.d0   ! adjust after a first test run

  allocate(t_field(nx,0:ny+1))
  allocate(phi_field(nx,0:ny+1))
  allocate(cond_phi_sum(n_bins), cond_count(n_bins))

  cond_phi_sum = 0.d0
  cond_count   = 0.d0
  bin_width    = gradT_max / dble(n_bins)
  n_snaps      = 0

  do snap = t_start, t_end, t_step
    write(fname_t,   '("../../src/output/t_",I8.8,".dat")')   snap
    write(fname_phi, '("../../src/output/phi_",I8.8,".dat")') snap

    open(unit=10, file=trim(fname_t),   form='unformatted', access='stream', status='old')
    read(10) t_field
    close(10)

    open(unit=11, file=trim(fname_phi), form='unformatted', access='stream', status='old')
    read(11) phi_field
    close(11)

    do j = 1, ny
      do i = 1, nx

        ! -----------------------------------------------------------
        ! x-gradient: central difference, periodic in x
        ! -----------------------------------------------------------
        if (i == 1) then
          dTdx = (t_field(2,j) - t_field(nx,j)) / (2.d0*dx)
        else if (i == nx) then
          dTdx = (t_field(1,j) - t_field(nx-1,j)) / (2.d0*dx)
        else
          dTdx = (t_field(i+1,j) - t_field(i-1,j)) / (2.d0*dx)
        end if

        ! -----------------------------------------------------------
        ! y-gradient: central difference, T is cell-centered
        ! ghost cells (j=0, j=ny+1) should be set by BCs
        ! For isothermal walls: T(i,0) = t_hot, T(i,ny+1) = t_cold
        ! (or vice versa depending on your convention)
        ! -----------------------------------------------------------
        dTdy = (t_field(i,j+1) - t_field(i,j-1)) / (2.d0*dy)

        gradT_mag = dsqrt(dTdx**2 + dTdy**2)

        ! Keep only points inside the chosen range
        if (gradT_mag >= 0.d0 .and. gradT_mag < gradT_max) then
          ibin = int(gradT_mag / bin_width) + 1
          ibin = max(1, min(n_bins, ibin))
          cond_phi_sum(ibin) = cond_phi_sum(ibin) + phi_field(i,j)
          cond_count(ibin)   = cond_count(ibin)   + 1.d0
        end if

      end do
    end do

    n_snaps = n_snaps + 1
    write(*,'("Processed snapshot ",I8)') snap
  end do

  open(unit=30, file='cond_avg_phi_gradT.dat', status='replace')
  write(30,'(A)') '# gradT_bin_center   <phi|gradT>   count'
  do ibin = 1, n_bins
    gradT_bin_center = (dble(ibin) - 0.5d0) * bin_width
    if (cond_count(ibin) > 0.d0) then
      phi_cond = cond_phi_sum(ibin) / cond_count(ibin)
    else
      phi_cond = 0.d0
    end if
    write(30,'(3(ES20.10,2X))') gradT_bin_center, phi_cond, cond_count(ibin)
  end do
  close(30)

  write(*,'("Written cond_avg_phi_gradT.dat (",I6," snapshots)")') n_snaps

  deallocate(t_field, phi_field)
  deallocate(cond_phi_sum, cond_count)

end program cond_avg_phi_gradT