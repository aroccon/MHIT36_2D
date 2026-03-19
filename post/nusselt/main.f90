program nusselt_rb
  implicit none

  integer,  parameter :: nx = 2048, ny = 2048
  real(8),  parameter :: lx = 2048.d0, ly = 2048.d0
  real(8),  parameter :: kappa = 0.166d0
  real(8),  parameter :: dy = ly / ny

  integer,  parameter :: t_start = 4000000
  integer,  parameter :: t_end   = 10000000
  integer,  parameter :: t_step  = 4000

  real(8), allocatable :: t_field(:,:), v_field(:,:), phi_field(:,:)

  ! total
  real(8), allocatable :: tmean(:), nuc(:), nut(:), nu_total(:)
  ! phase 1: continuous (phi=0)
  real(8), allocatable :: tmean1(:), nuc1(:), nut1(:), nu_total1(:)
  ! phase 2: drops (phi=1)
  real(8), allocatable :: tmean2(:), nuc2(:), nut2(:), nu_total2(:)
  ! drop fraction per plane
  real(8), allocatable :: drop_frac(:)

  integer :: i, j, snap
  real(8) :: v_center, phi_sharp
  character(len=64) :: fname_t, fname_v, fname_phi, fname_out

  allocate( t_field  (nx, 0:ny+1) )
  allocate( v_field  (nx,   ny+1) )
  allocate( phi_field(nx, 0:ny+1) )

  allocate( tmean(0:ny+1), nuc(0:ny+1), nut(0:ny+1), nu_total(0:ny+1) )
  allocate( tmean1(0:ny+1), nuc1(0:ny+1), nut1(0:ny+1), nu_total1(0:ny+1) )
  allocate( tmean2(0:ny+1), nuc2(0:ny+1), nut2(0:ny+1), nu_total2(0:ny+1) )
  allocate( drop_frac(0:ny+1) )

  do snap = t_start, t_end, t_step

    ! --- Read fields
    write(fname_t,   '("t_",   I8.8, ".dat")') snap
    write(fname_v,   '("v_",   I8.8, ".dat")') snap
    write(fname_phi, '("phi_", I8.8, ".dat")') snap

    open(unit=10, file='../../src/output/' // trim(fname_t),   form='unformatted', access='stream', status='old')
    read(10) t_field;   close(10)
    open(unit=11, file='../../src/output/' // trim(fname_v),   form='unformatted', access='stream', status='old')
    read(11) v_field;   close(11)
    open(unit=12, file='../../src/output/' // trim(fname_phi), form='unformatted', access='stream', status='old')
    read(12) phi_field; close(12)

    ! --- Reset all arrays
    tmean=0.d0;  nuc=0.d0;  nut=0.d0;  nu_total=0.d0
    tmean1=0.d0; nuc1=0.d0; nut1=0.d0; nu_total1=0.d0
    tmean2=0.d0; nuc2=0.d0; nut2=0.d0; nu_total2=0.d0
    drop_frac=0.d0

    ! --- Plane averages
    do j = 1, ny
      do i = 1, nx
        ! sharp indicator: 1 if continuous (phi<0.5), 0 if drop (phi>0.5)
        phi_sharp = merge(1.d0, 0.d0, phi_field(i,j) < 0.5d0)

        tmean(j)  = tmean(j)  + t_field(i,j)
        tmean1(j) = tmean1(j) + t_field(i,j) * phi_sharp          ! continuous
        tmean2(j) = tmean2(j) + t_field(i,j) * (1.d0 - phi_sharp) ! drops
        drop_frac(j) = drop_frac(j) + (1.d0 - phi_sharp)
      end do

      tmean(j)     = tmean(j)     / nx
      drop_frac(j) = drop_frac(j) / nx   ! fraction of plane occupied by drops

      ! phase-averaged T: <indicator*T> / <indicator>
      ! continuous phase
      if ((1.d0 - drop_frac(j)) > 1.d-10) then
        tmean1(j) = tmean1(j) / nx / (1.d0 - drop_frac(j))
      else
        tmean1(j) = 0.d0
      end if
      ! drops
      if (drop_frac(j) > 1.d-10) then
        tmean2(j) = tmean2(j) / nx / drop_frac(j)
      else
        tmean2(j) = 0.d0
      end if
    end do

    ! --- Conductive Nu (centered difference, weighted by phase fraction)
    do j = 1, ny
      nuc(j)  = -0.5d0*(tmean(j+1)  - tmean(j-1) ) / dy * ly
      ! weight by phase fraction so nuc1 + nuc2 = nuc
      nuc1(j) = -0.5d0*(tmean1(j+1) - tmean1(j-1)) / dy * ly * (1.d0 - drop_frac(j))
      nuc2(j) = -0.5d0*(tmean2(j+1) - tmean2(j-1)) / dy * ly * drop_frac(j)
    end do

    ! --- Convective Nu (sharp phase discrimination)
    do j = 1, ny
      do i = 1, nx
        v_center  = 0.5d0 * (v_field(i,j) + v_field(i,j+1))
        phi_sharp = merge(1.d0, 0.d0, phi_field(i,j) < 0.5d0)
        nut(j)    = nut(j)  + t_field(i,j) * v_center
        nut1(j)   = nut1(j) + t_field(i,j) * v_center * phi_sharp          ! continuous
        nut2(j)   = nut2(j) + t_field(i,j) * v_center * (1.d0 - phi_sharp) ! drops
      end do
      nut(j)  = nut(j)  / nx * ly / kappa
      nut1(j) = nut1(j) / nx * ly / kappa
      nut2(j) = nut2(j) / nx * ly / kappa
    end do

    ! --- Total Nu per phase
    do j = 1, ny
      nu_total(j)  = nuc(j)  + nut(j)
      nu_total1(j) = nuc1(j) + nut1(j)
      nu_total2(j) = nuc2(j) + nut2(j)
    end do

    ! --- Write output ASCII file
    write(fname_out, '("nusselt_", I8.8, ".dat")') snap
    open(unit=20, file=trim(fname_out), status='replace')
    write(20,'(A)') '# j  drop_frac  nuc  nut  nu_total  nuc1_cont  nut1_cont  nutotal1_cont  nuc2_drop  nut2_drop  nutotal2_drop'
    do j = 1, ny
      write(20,'(I6, 10(2X,ES15.7))') j, drop_frac(j), &
            nuc(j),  nut(j),  nu_total(j),  &
            nuc1(j), nut1(j), nu_total1(j), &
            nuc2(j), nut2(j), nu_total2(j)
    end do
    close(20)

    write(*,'("Processed snapshot ",I8)') snap

  end do

  deallocate(t_field, v_field, phi_field)
  deallocate(tmean, nuc, nut, nu_total)
  deallocate(tmean1, nuc1, nut1, nu_total1)
  deallocate(tmean2, nuc2, nut2, nu_total2)
  deallocate(drop_frac)

end program nusselt_rb
