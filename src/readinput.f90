subroutine readinput
use velocity
use phase
use param
implicit none

open(unit=55,file='input.inp',form='formatted',status='old')
!Time step parameters
read(55,*) restart
read(55,*) tstart
read(55,*) tfin
read(55,*) dump
!Flow parameters
read(55,*) dt
read(55,*) rho
read(55,*) mu
! temperature parameters
read(55,*) alphag
read(55,*) pr
read(55,*) lx
read(55,*) ly
read(55,*) utop
read(55,*) ubot
! phase-field parameters
read(55,*) icphi
read(55,*) radius
read(55,*) sigma
read(55,*) epsr
read(55,*) xoffset
read(55,*) yoffset
! near-contact repulsive force (Liu et al., PoF 37, 092123, 2025)
read(55,*) ahamaker
read(55,*) hcr
read(55,*) hminr
read(55,*) hmaxr



! compute pre-defined constant 
twopi=2.d0*pi 
difftemp=mu/pr
!difftemp=sqrt(pr/ra)! sqrt(Ra) with Ra=2e6 
!mu=sqrt(1.d0/ra) ! sqrt(1/Ra) with Ra=2e6 
dx = lx/nx
dy = ly/ny
dxi=1.d0/dx
ddxi=1.d0/dx/dx
dyi=1.d0/dy
ddyi=1.d0/dy/dy
rhoi=1.d0/rho
eps=epsr*max(dx,dy)
epsi=1.d0/eps
hc=hcr*max(dx,dy)
hmin=hminr*max(dx,dy)
hmax=hmaxr*max(dx,dy)
enum=1.e-16
tbot=0.5d0
ttop=-0.5d0


write(*,*) "------------------------------------------------------"
write(*,*) "@@@@@@@@@@  @@@@@@@  @@@@@@@   @@@@@@@ @@@@@@    @@@@@"
write(*,*) "@@! @@! @@! @@!  @@@ @@!  @@@ !@@          @@! @@!@"    
write(*,*) "@!! !!@ @!@ @!@!!@!  @!@!@!@  !@!       @!!!:  @!@!@!@"
write(*,*) "!!:     !!: !!: :!!  !!:  !!! :!!          !!: !!:  !!!"
write(*,*) " :      :    :   : : :: : ::   :: :: : ::: ::   : : ::"
write(*,*) "------------------------------------------------------"
write(*,*) "Grid:    ", nx, 'x', ny
write(*,*) "Tfin     ", tfin
write(*,*) "Dump        ", dump
write(*,*) "Density      ", rho
write(*,*) "Viscosity      ", mu
write(*,*) "Alphag        ", alphag
write(*,*) "Prandtl        ", pr
write(*,*) "Difftemp (from Pr)", difftemp
write(*,*) "Radius          ", radius
if (icphi .eq. 2) write(*,*) "Drop offsets x,y", xoffset, yoffset
write(*,*) "Sigma          ", sigma
write(*,*) "Eps             ", eps
write(*,*) "Hamaker A_H     ", ahamaker
write(*,*) "hc, hmin, hmax  ", hc, hmin, hmax
write(*,*) 'Lx             ', lx
write(*,*) 'Ly             ', ly
write(*,*) 'U top wall     ', utop
write(*,*) 'U bottom wall  ', ubot
write(*,*) 'Dx              ', dx
write(*,*) 'Dy              ', dy
end subroutine


