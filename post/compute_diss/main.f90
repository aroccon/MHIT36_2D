program dissipation_rb
  implicit none

  integer,  parameter :: nx = 2048, ny = 2048
  real(8),  parameter :: lx = 2048.d0, ly = 2048.d0
  real(8),  parameter :: dx = lx / dble(nx)
  real(8),  parameter :: dy = ly / dble(ny)
  real(8),  parameter :: nu_visc = 0.166d0  
  integer,  parameter :: t_start = 4000000
  integer,  parameter :: t_end   = 8000000
  integer,  parameter :: t_step  = 4000
  ! u: staggered in x -> u(0:nx, ny)  face at i+1/2
  ! v: staggered in y -> v(nx, 0:ny)  face at j+1/2
  real(8), allocatable :: u_field(:,:)
  real(8), allocatable :: v_field(:,:)
  integer :: i, j, snap, n_snaps, im, ip, jm, jp
  real(8) :: dudx, dudy, dvdx, dvdy
  real(8) :: eps_u, eps_u_run

  character(len=64) :: fname_u, fname_v

  allocate(u_field(nx,0:ny+1))
  allocate(v_field(nx,ny+1))

  n_snaps   = 0
  eps_u_run = 0.d0

  open(unit=30, file='dissipation_timeseries.dat', status='replace')
  write(30,'(A)') '# snap   eps_u   eps_u_running_avg'

  do snap = t_start, t_end, t_step

    ! read fields
    write(fname_u, '("u_", I8.8, ".dat")') snap
    write(fname_v, '("v_", I8.8, ".dat")') snap
    open(unit=11, file='../../src/output/' // trim(fname_u), form='unformatted', access='stream', status='old')
    read(11) u_field; close(11)
    open(unit=12, file='../../src/output/' // trim(fname_v), form='unformatted', access='stream', status='old')
    read(12) v_field; close(12)

    eps_u = 0.d0

    ! compute dissipation at cell center
    do j = 1, ny
      do i = 1, nx
        ip = i + 1
        im = i - 1
        jp = j + 1
        jm = j - 1
        ! Periodic wrapping in x
        if (im < 1)  im = nx
        if (ip > nx) ip = 1
        ! dudx: direct from staggered u-faces
        dudx = (u_field(ip,j) - u_field(i,j)) / dx
        ! dvdy: direct from staggered v-faces
        dvdy = (v_field(i,jp) - v_field(i,j)) / dy
        ! dudy: interpolate u to cell corners then difference
        dudy = 0.5d0 * ( (u_field(i,jp) + u_field(ip,jp)) &
                       - (u_field(i,jm) + u_field(ip,jm)) )/(2.d0*dy)
        ! dvdx: interpolate v to cell corners then difference
        dvdx = 0.5d0 * ( (v_field(ip,j) + v_field(ip,jp)) &
                       - (v_field(im,j) + v_field(im,jp)) )/(2.d0*dx)
        eps_u = eps_u + nu_visc * ( 2.d0*(dudx**2 + dvdy**2)  + (dudy + dvdx)**2 )
      end do
    end do

    ! --- Normalize by number of cells
    eps_u = eps_u / dble(nx*ny)

    ! --- Update running average
    n_snaps   = n_snaps + 1
    eps_u_run = eps_u_run + (eps_u - eps_u_run) / dble(n_snaps)

    ! --- Write to file
    write(30, '(I10, 2(2X,ES15.7))') snap, eps_u, eps_u_run
    write(*,  '("Snapshot ",I8," eps_u=",ES15.7," <eps_u>=",ES15.7)') &
              snap, eps_u, eps_u_run

  end do

  close(30)
  write(*,'("Done. Written dissipation_timeseries.dat (",I0," snapshots)")') n_snaps

  deallocate(u_field, v_field)

end program dissipation_rb