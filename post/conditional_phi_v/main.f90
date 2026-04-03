program cond_avg_phi_V
  implicit none
  integer,  parameter :: nx = 2048, ny = 2048
  integer,  parameter :: t_start = 2000000
  integer,  parameter :: t_end   = 10000000
  integer,  parameter :: t_step  = 4000
  integer,  parameter :: n_vbins = 100
  real(8),  parameter :: v_min   = -0.5d0   ! adjust to your expected velocity range
  real(8),  parameter :: v_max   =  0.5d0   ! adjust to your expected velocity range

  real(8), allocatable :: v_field(:,:), phi_field(:,:)
  real(8), allocatable :: cond_phi_sum(:), cond_count(:)

  integer :: i, j, snap, ibin, n_snaps
  real(8) :: v_center, bin_width, v_bin_center, phi_cond
  character(len=128) :: fname_v, fname_phi

  allocate(v_field(nx,1:ny+1))
  allocate(phi_field(nx,0:ny+1))
  allocate(cond_phi_sum(n_vbins), cond_count(n_vbins))

  cond_phi_sum = 0.d0
  cond_count   = 0.d0
  bin_width    = (v_max - v_min) / dble(n_vbins)
  n_snaps      = 0

  do snap = t_start, t_end, t_step
    write(fname_v,   '("../../src/output/v_",I8.8,".dat")')   snap
    write(fname_phi, '("../../src/output/phi_",I8.8,".dat")') snap

    open(unit=10, file=trim(fname_v), form='unformatted', access='stream', status='old')
    read(10) v_field
    close(10)

    open(unit=11, file=trim(fname_phi), form='unformatted', access='stream', status='old')
    read(11) phi_field
    close(11)

    do j = 1, ny
      do i = 1, nx
        ! Interpolate staggered v to cell center: average v(i,j) and v(i,j+1)
        v_center = 0.5d0 * (v_field(i,j) + v_field(i,j+1))

        ! Keep only points inside the chosen velocity range
        if (v_center >= v_min .and. v_center < v_max) then
          ibin = int((v_center - v_min) / bin_width) + 1
          ! Safety clamp
          ibin = max(1, min(n_vbins, ibin))
          cond_phi_sum(ibin) = cond_phi_sum(ibin) + phi_field(i,j)
          cond_count(ibin)   = cond_count(ibin)   + 1.d0
        end if
      end do
    end do

    n_snaps = n_snaps + 1
    write(*,'("Processed snapshot ",I8)') snap
  end do

  open(unit=30, file='cond_avg_phi_V.dat', status='replace')
  write(30,'(A)') '# V_bin_center   <phi|V>   count'
  do ibin = 1, n_vbins
    v_bin_center = v_min + (dble(ibin) - 0.5d0) * bin_width
    if (cond_count(ibin) > 0.d0) then
      phi_cond = cond_phi_sum(ibin) / cond_count(ibin)
    else
      phi_cond = 0.d0
    end if
    write(30,'(3(ES20.10,2X))') v_bin_center, phi_cond, cond_count(ibin)
  end do
  close(30)

  write(*,'("Written cond_avg_phi_V.dat (",I6," snapshots)")') n_snaps

  deallocate(v_field, phi_field)
  deallocate(cond_phi_sum, cond_count)

end program cond_avg_phi_V