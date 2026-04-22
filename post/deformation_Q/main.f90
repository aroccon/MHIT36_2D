program drop_count
  implicit none

  integer, parameter :: nx = 2048, ny = 2048
  integer, parameter :: t_start = 4000000
  integer, parameter :: t_end   = 8000000
  integer, parameter :: t_step  = 4000

  real(8), parameter :: cell_area = 1.d0
  real(8), parameter :: lx = 2048.d0, ly = 2048.d0
  real(8), parameter :: ddx = lx / dble(nx)
  real(8), parameter :: ddy = ly / dble(ny)

  ! 8-connected neighbour offsets
  integer, parameter :: n_nb = 8
  integer, parameter :: dx(8) = (/ -1, 0, 1, -1, 1, -1, 0, 1 /)
  integer, parameter :: dy(8) = (/ -1,-1,-1,  0, 0,  1, 1, 1 /)

  ! --- Field storage
  ! phi : cell-centered, with y-ghost cells (j=0 and j=ny+1)
  ! u   : staggered on LEFT face of cell (i,j); shape (nx, 0:ny+1),
  !       j=0 and j=ny+1 are ghost rows for no-slip at walls.
  ! v   : staggered on BOTTOM face of cell (i,j); shape (nx, 1:ny+1),
  !       v(:,1) is on the bottom wall and v(:,ny+1) on the top wall
  !       (no ghost cells needed, walls are real face values).
  real(8), allocatable :: phi_field(:,:)
  real(8), allocatable :: u_field(:,:), v_field(:,:)
  real(8), allocatable :: uc(:,:), vc(:,:)    ! cell-centered velocities
  real(8), allocatable :: Q_field(:,:)
  logical, allocatable :: mask(:,:)
  logical, allocatable :: visited(:,:)

  ! stack for flood fill
  integer, allocatable :: stack_x(:), stack_y(:), stack_ux(:)
  integer :: stack_top

  ! per-drop arrays
  integer, parameter   :: max_drops = 100000
  integer, allocatable :: drop_area(:)
  integer, allocatable :: drop_xmin(:), drop_xmax(:)
  integer, allocatable :: drop_ymin(:), drop_ymax(:)
  integer, allocatable :: drop_uxmin(:), drop_uxmax(:)
  real(8), allocatable :: drop_Qsum(:)
  logical, allocatable :: drop_wraps(:)
  real(8), allocatable :: drop_xsum(:), drop_ysum(:)
  real(8), allocatable :: drop_Ixx(:), drop_Iyy(:), drop_Ixy(:)

  ! local variables
  integer :: i, j, snap, n_drops, d
  integer :: ci, cj, ni, nj
  integer :: cux, nux
  integer :: im1, ip1
  real(8) :: xc, yc, Ixx, Iyy, Ixy
  real(8) :: trace, det, disc, lambda1, lambda2, deformation
  real(8) :: dudx, dudy, dvdx, dvdy
  real(8) :: Sxx, Syy, Sxy, Oxy, S2, O2
  character(len=64) :: fname_phi, fname_u, fname_v, fname_out

  allocate( phi_field(nx, 0:ny+1) )
  allocate( u_field  (nx, 0:ny+1) )    ! left-face u, with y-ghosts
  allocate( v_field  (nx, 1:ny+1) )    ! bottom-face v, walls included
  allocate( uc       (nx, ny) )        ! cell-centered u
  allocate( vc       (nx, ny) )        ! cell-centered v
  allocate( Q_field  (nx, ny) )
  allocate( mask(nx, ny), visited(nx, ny) )
  allocate( stack_x(nx*ny), stack_y(nx*ny), stack_ux(nx*ny) )
  allocate( drop_area(max_drops) )
  allocate( drop_xmin(max_drops), drop_xmax(max_drops) )
  allocate( drop_ymin(max_drops), drop_ymax(max_drops) )
  allocate( drop_uxmin(max_drops), drop_uxmax(max_drops) )
  allocate( drop_Qsum(max_drops) )
  allocate( drop_wraps(max_drops) )
  allocate( drop_xsum(max_drops), drop_ysum(max_drops) )
  allocate( drop_Ixx(max_drops), drop_Iyy(max_drops), drop_Ixy(max_drops) )

  do snap = t_start, t_end, t_step

    ! --- Read fields
    write(fname_phi, '("phi_", I8.8, ".dat")') snap
    write(fname_u,   '("u_",   I8.8, ".dat")') snap
    write(fname_v,   '("v_",   I8.8, ".dat")') snap

    open(unit=10, file='../../src/output/' // trim(fname_phi), &
         form='unformatted', access='stream', status='old')
    read(10) phi_field
    close(10)

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

    ! --- Build binary mask
    do j = 1, ny
      do i = 1, nx
        mask(i,j) = ( phi_field(i,j) > 0.5d0 )
      end do
    end do

    visited = .false.
    n_drops = 0

    ! --- Scan and flood fill
    do j = 1, ny
      do i = 1, nx

        if ( mask(i,j) .and. .not. visited(i,j) ) then

          n_drops = n_drops + 1
          if (n_drops > max_drops) then
            write(*,*) 'ERROR: max_drops exceeded at snapshot ', snap
            stop
          end if

          ! initialise drop properties
          drop_area(n_drops)  = 0
          drop_Qsum(n_drops)  = 0.d0
          drop_xmin(n_drops)  = i
          drop_xmax(n_drops)  = i
          drop_ymin(n_drops)  = j
          drop_ymax(n_drops)  = j
          drop_uxmin(n_drops) = i
          drop_uxmax(n_drops) = i
          drop_wraps(n_drops) = .false.
          drop_xsum(n_drops)  = 0.d0
          drop_ysum(n_drops)  = 0.d0
          drop_Ixx(n_drops)   = 0.d0
          drop_Iyy(n_drops)   = 0.d0
          drop_Ixy(n_drops)   = 0.d0

          ! push seed
          stack_top    = 1
          stack_x(1)   = i
          stack_y(1)   = j
          stack_ux(1)  = i
          visited(i,j) = .true.

          do while (stack_top > 0)

            ! pop
            ci  = stack_x(stack_top)
            cj  = stack_y(stack_top)
            cux = stack_ux(stack_top)
            stack_top = stack_top - 1

            drop_area(n_drops) = drop_area(n_drops) + 1

            ! --- Flow topology parameter at this cell
            drop_Qsum(n_drops) = drop_Qsum(n_drops) + Q_field(ci,cj)

            ! accumulate position sums for centroid (unwrapped x)
            drop_xsum(n_drops) = drop_xsum(n_drops) + dble(cux)
            drop_ysum(n_drops) = drop_ysum(n_drops) + dble(cj)

            drop_uxmin(n_drops) = min(drop_uxmin(n_drops), cux)
            drop_uxmax(n_drops) = max(drop_uxmax(n_drops), cux)
            drop_ymin(n_drops)  = min(drop_ymin(n_drops),  cj)
            drop_ymax(n_drops)  = max(drop_ymax(n_drops),  cj)

            ! visit 8 neighbours
            do d = 1, n_nb
              ni  = mod(ci + dx(d) - 1 + nx, nx) + 1
              nj  = cj + dy(d)
              nux = cux + dx(d)

              if (nj < 1 .or. nj > ny) cycle

              if ( mask(ni,nj) .and. .not. visited(ni,nj) ) then
                visited(ni,nj) = .true.
                stack_top = stack_top + 1
                stack_x(stack_top)  = ni
                stack_y(stack_top)  = nj
                stack_ux(stack_top) = nux
              end if
            end do

          end do ! flood fill

          ! --- Convert unwrapped x-bounds back to periodic indices
          drop_xmin(n_drops) = mod(drop_uxmin(n_drops) - 1 + nx, nx) + 1
          drop_xmax(n_drops) = mod(drop_uxmax(n_drops) - 1 + nx, nx) + 1
          drop_wraps(n_drops) = (drop_uxmin(n_drops) < 1) .or. &
                                (drop_uxmax(n_drops) > nx)

          ! --- Compute centroid using unwrapped x
          xc = drop_xsum(n_drops) / dble(drop_area(n_drops))
          yc = drop_ysum(n_drops) / dble(drop_area(n_drops))

          ! --- Second pass over bounding box: compute inertia tensor
          do cj = drop_ymin(n_drops), drop_ymax(n_drops)
            do cux = drop_uxmin(n_drops), drop_uxmax(n_drops)
              ci = mod(cux - 1 + nx, nx) + 1
              if (mask(ci,cj)) then
                drop_Ixx(n_drops) = drop_Ixx(n_drops) + (dble(cj)  - yc)**2
                drop_Iyy(n_drops) = drop_Iyy(n_drops) + (dble(cux) - xc)**2
                drop_Ixy(n_drops) = drop_Ixy(n_drops) + (dble(cux) - xc) &
                                                       * (dble(cj)  - yc)
              end if
            end do
          end do

          ! --- Normalize by area
          Ixx = drop_Ixx(n_drops) / dble(drop_area(n_drops))
          Iyy = drop_Iyy(n_drops) / dble(drop_area(n_drops))
          Ixy = drop_Ixy(n_drops) / dble(drop_area(n_drops))

          ! --- Eigenvalues of 2x2 symmetric inertia tensor
          trace   = Ixx + Iyy
          det     = Ixx*Iyy - Ixy**2
          disc    = dsqrt(max(0.d0, 0.25d0*trace**2 - det))
          lambda1 = 0.5d0*trace + disc
          lambda2 = 0.5d0*trace - disc

          ! --- Deformation parameter
          if (lambda1 + lambda2 > 0.d0) then
            deformation = (lambda1 - lambda2) / (lambda1 + lambda2)
          else
            deformation = 0.d0
          end if

          ! --- Store in arrays for output
          drop_Ixx(n_drops) = lambda1
          drop_Iyy(n_drops) = lambda2
          drop_Ixy(n_drops) = deformation

        end if

      end do
    end do

    ! --- Write output
    write(fname_out, '("drops_", I8.8, ".dat")') snap
    open(unit=20, file=trim(fname_out), status='replace')
    write(20,'(A)') '% drop_id  area  x_min  x_max  y_min  y_max  ' // &
                    'Q_mean  lambda1  lambda2  deformation'
    do d = 1, n_drops
      write(20,'(I6, 2X, I8, 4(2X,I6), 4(2X,ES15.7))') d,            &
        drop_area(d),                                                   &
        drop_xmin(d), drop_xmax(d),                                    &
        drop_ymin(d), drop_ymax(d),                                    &
        drop_Qsum(d) / dble(drop_area(d)),                             &   ! Q_mean
        drop_Ixx(d),                                                   &   ! lambda1
        drop_Iyy(d),                                                   &   ! lambda2
        drop_Ixy(d)                                                        ! deformation
    end do
    close(20)

    write(*,'("Snapshot ",I8,": ",I6," drops found")') snap, n_drops

  end do

  deallocate(phi_field, u_field, v_field, uc, vc, Q_field)
  deallocate(mask, visited)
  deallocate(stack_x, stack_y, stack_ux)
  deallocate(drop_area, drop_xmin, drop_xmax, drop_ymin, drop_ymax)
  deallocate(drop_uxmin, drop_uxmax)
  deallocate(drop_Qsum, drop_wraps)
  deallocate(drop_xsum, drop_ysum)
  deallocate(drop_Ixx, drop_Iyy, drop_Ixy)

end program drop_count