program cond_avg_phi_T
  implicit none

  integer,  parameter :: nx = 2048, ny = 2048
  integer,  parameter :: t_start = 2000000
  integer,  parameter :: t_end   = 10000000
  integer,  parameter :: t_step  = 4000

  integer,  parameter :: n_tbins = 100
  real(8),  parameter :: t_hot   =  0.5d0
  real(8),  parameter :: t_cold  = -0.5d0

  real(8), allocatable :: t_field(:,:), phi_field(:,:)
  real(8), allocatable :: cond_phi_sum(:), cond_count(:)

  integer :: i, j, snap, ibin, n_snaps
  real(8) :: temp, bin_width, t_bin_center, phi_cond
  character(len=128) :: fname_t, fname_phi

  allocate(t_field(nx,0:ny+1))
  allocate(phi_field(nx,0:ny+1))
  allocate(cond_phi_sum(n_tbins), cond_count(n_tbins))

  cond_phi_sum = 0.d0
  cond_count   = 0.d0
  bin_width    = (t_hot - t_cold) / dble(n_tbins)
  n_snaps      = 0

  do snap = t_start, t_end, t_step

    write(fname_t,   '("../../src/output/t_",I8.8,".dat")')   snap
    write(fname_phi, '("../../src/output/phi_",I8.8,".dat")') snap

    open(unit=10, file=trim(fname_t), form='unformatted', access='stream', status='old')
    read(10) t_field
    close(10)

    open(unit=11, file=trim(fname_phi), form='unformatted', access='stream', status='old')
    read(11) phi_field
    close(11)

    do j = 1, ny
      do i = 1, nx
        temp = t_field(i,j)
        ! Keep only points inside the chosen temperature range
        if (temp >= t_cold .and. temp < t_hot) then
          ibin = int((temp - t_cold) / bin_width) + 1

          ! Safety clamp
          ibin = max(1, min(n_tbins, ibin))

          cond_phi_sum(ibin) = cond_phi_sum(ibin) + phi_field(i,j)
          cond_count(ibin)   = cond_count(ibin)   + 1.d0
        end if

      end do
    end do

    n_snaps = n_snaps + 1
    write(*,'("Processed snapshot ",I8)') snap

  end do

  open(unit=30, file='cond_avg_phi_T.dat', status='replace')
  write(30,'(A)') '# T_bin_center   <phi|T>   count'

  do ibin = 1, n_tbins
    t_bin_center = t_cold + (dble(ibin) - 0.5d0) * bin_width

    if (cond_count(ibin) > 0.d0) then
      phi_cond = cond_phi_sum(ibin) / cond_count(ibin)
    else
      phi_cond = 0.d0
    end if

    write(30,'(3(ES20.10,2X))') t_bin_center, phi_cond, cond_count(ibin)
  end do

  close(30)

  write(*,'("Written cond_avg_phi_T.dat (",I6," snapshots)")') n_snaps

  deallocate(t_field, phi_field)
  deallocate(cond_phi_sum, cond_count)

end program cond_avg_phi_T