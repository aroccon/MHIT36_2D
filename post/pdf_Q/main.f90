program q_pdf
  implicit none

  integer, parameter :: nx = 2048, ny = 2048
  integer, parameter :: t_start = 4000000
  integer, parameter :: t_end   = 8000000
  integer, parameter :: t_step  = 4000

  real(8), parameter :: lx = 2048.d0, ly = 2048.d0
  real(8), parameter :: ddx = lx / dble(nx)
  real(8), parameter :: ddy = ly / dble(ny)

  ! --- PDF binning
  integer, parameter :: n_bins = 100
  real(8), parameter :: q_min = -1.d0
  real(8), parameter :: q_max =  1.d0
  real(8), parameter :: db    = (q_max - q_min) / dble(n_bins)

  ! --- Field storage
  ! u : staggered on LEFT face of cell (i,j); shape (nx, 0:ny+1),
  !     j=0 and j=ny+1 are ghost rows for no-slip at walls.
  ! v : staggered on BOTTOM face of cell (i,j); shape (nx, 1:ny+1),
  !     v(:,1) is on the bottom wall and v(:,ny+1) on the top wall.
  real(8), allocatable :: u_field(:,:), v_field(:,:)
  real(8), allocatable :: uc(:,:), vc(:,:)    ! cell-centered velocities
  real(8), allocatable :: Q_field(:,:)

  ! PDF accumulators (global over all snapshots)
  integer(8), allocatable :: hist(:)
  integer(8) :: total_count
  real(8), allocatable :: pdf(:), bin_center(:)

  ! local variables
  integer :: i, j, snap, b
  integer :: im1, ip1
  real(8) :: dudx, dudy, dvdx, dvdy
  real(8) :: Sxx, Syy, Sxy, Oxy, S2, O2, Qval
  character(len=64) :: fname_u, fname_v, fname_out

  allocate( u_field (nx, 0:ny+1) )    ! left-face u, with y-ghosts
  allocate( v_field (nx, 1:ny+1) )    ! bottom-face v, walls included
  allocate( uc      (nx, ny) )        ! cell-centered u
  allocate( vc      (nx, ny) )        ! cell-centered v
  allocate( Q_field (nx, ny) )
  allocate( hist(n_bins) )
  allocate( pdf(n_bins), bin_center(n_bins) )

  hist = 0_8
  total_count = 0_8

  do b = 1, n_bins
    bin_center(b) = q_min + (dble(b) - 0.5d0)*db
  end do

  do snap = t_start, t_end, t_step

    ! --- Read fields
    write(fname_u, '("u_", I8.8, ".dat")') snap
    write(fname_v, '("v_", I8.8, ".dat")') snap

    open(unit=12, file='../../src/output/' // trim(fname_u), &
         form='unformatted', access='stream', status='old')
    read(12) u_field
    close(12)

    open(unit=13, file='../../src/output/' // trim(fname_v), &
         form='unformatted', access='stream', status='old')
    read(13) v_field
    close(13)

    ! --- Interpolate staggered velocities to cell centers
    !     u(i,j) on LEFT face of cell (i,j):   uc(i,j) = 0.5*(u(i,j) + u(i+1,j))
    !     v(i,j) on BOTTOM face of cell (i,j): vc(i,j) = 0.5*(v(i,j) + v(i,j+1))
    !     For u, i+1 wraps periodically in x.
    do j = 1, ny
      do i = 1, nx
        ip1 = mod(i, nx) + 1
        uc(i,j) = 0.5d0 * ( u_field(i,  j) + u_field(ip1, j) )
        vc(i,j) = 0.5d0 * ( v_field(i,  j) + v_field(i,   j+1) )
      end do
    end do

    ! --- Flow topology parameter Q at cell centers
    !     Q = (||S||^2 - ||Omega||^2) / (||S||^2 + ||Omega||^2)
    !     Q -> -1 : pure rotation
    !     Q ->  0 : pure shear
    !     Q -> +1 : pure extension
    !
    !     Gradients: periodic central differences in x,
    !     one-sided at j=1 and j=ny, central elsewhere.
    do j = 1, ny
      do i = 1, nx
        im1 = mod(i-2+nx, nx) + 1
        ip1 = mod(i,      nx) + 1

        dudx = (uc(ip1,j) - uc(im1,j)) / (2.d0*ddx)
        dvdx = (vc(ip1,j) - vc(im1,j)) / (2.d0*ddx)

        if (j == 1) then
          dudy = (uc(i,j+1) - uc(i,j)) / ddy
          dvdy = (vc(i,j+1) - vc(i,j)) / ddy
        else if (j == ny) then
          dudy = (uc(i,j) - uc(i,j-1)) / ddy
          dvdy = (vc(i,j) - vc(i,j-1)) / ddy
        else
          dudy = (uc(i,j+1) - uc(i,j-1)) / (2.d0*ddy)
          dvdy = (vc(i,j+1) - vc(i,j-1)) / (2.d0*ddy)
        end if

        Sxx = dudx
        Syy = dvdy
        Sxy = 0.5d0*(dudy + dvdx)
        Oxy = 0.5d0*(dudy - dvdx)

        S2 = Sxx**2 + Syy**2 + 2.d0*Sxy**2
        O2 = 2.d0*Oxy**2

        if (S2 + O2 > 1.d-30) then
          Q_field(i,j) = (S2 - O2) / (S2 + O2)
        else
          Q_field(i,j) = 0.d0
        end if
      end do
    end do

    ! --- Accumulate Q values into histogram
    do j = 1, ny
      do i = 1, nx
        Qval = Q_field(i,j)
        b = int( (Qval - q_min) / db ) + 1
        if (b < 1)      b = 1
        if (b > n_bins) b = n_bins
        hist(b)     = hist(b) + 1_8
        total_count = total_count + 1_8
      end do
    end do

    write(*,'("Snapshot ",I8," processed")') snap

  end do

  ! --- Normalize: PDF such that sum(pdf)*db = 1
  do b = 1, n_bins
    pdf(b) = dble(hist(b)) / ( dble(total_count) * db )
  end do

  ! --- Write output
  fname_out = 'Q_pdf.dat'
  open(unit=20, file=trim(fname_out), status='replace')
  write(20,'(A)') '% bin_center  pdf  count'
  do b = 1, n_bins
    write(20,'(ES15.7, 2X, ES15.7, 2X, I12)') bin_center(b), pdf(b), hist(b)
  end do
  close(20)

  write(*,'("Total samples: ",I0)') total_count
  write(*,'("PDF written to ",A)') trim(fname_out)

  deallocate(u_field, v_field, uc, vc, Q_field)
  deallocate(hist, pdf, bin_center)

end program q_pdf