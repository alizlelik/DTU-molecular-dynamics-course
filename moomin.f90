Module ParDataStructure

      double precision :: avogad=6.022045e+23       ! (# of atoms)/mol
      double precision :: boltz=1.380662e-23        ! J/K
	  
	  double precision :: rho, time_simulated, L, T_system, Lred
	  double precision :: p_system, deltar, pi !these are the new ones introduced by me
	  
	  
      double precision :: epsk,eps,epsred,sigma,xmass,mu,timeps
      double precision :: dt,dtsq
      double precision :: einit,ered,req,velinit
      double precision :: x,x0,xnew,xold,vx,acc
      double precision :: Epot,Ekin,hamilt

      integer :: maxstep,imd
	  
	  integer :: N !these are the new ones 
	  
	  logical :: isTailOn
	  
	  double precision, allocatable, dimension(:) :: h_t
	  integer, dimension(60) :: g_hist
	  character(3) :: ens, config
End Module ParDataStructure



Module AtomsDataStructure

	double precision, allocatable, dimension(:) :: coordx, coordy, coordz
	double precision, allocatable, dimension(:) :: coordx_old, coordy_old, coordz_old
	double precision, allocatable, dimension(:) :: velx, vely, velz
	double precision, allocatable, dimension(:) :: Fx,Fy,Fz
	double precision, allocatable, dimension(:,:) :: forcesx, forcesy, forcesz !stored only for pressure calc
	!double precision, allocatable, dimension(:) :: accx,accy,accz !old accelerations for velocity verlet
	

End Module AtomsDataStructure

module NVTDataStructure
	double precision :: Ms, taus, zeta, Tdes, lns, z, comp
	double precision, allocatable, dimension(:) :: z_t
end module

module accelerations !moved  here for velocity verlet debugging
	double precision, allocatable, dimension(:) :: accx,accy,accz

end module

subroutine gaussian(delta2,d,x1,x2) !generate the numbers with a Gaussian distribution
implicit none
double precision, intent(in) :: delta2 !read the deltasquared
integer, intent(in) :: d !read how many numbers we want
double precision, dimension(d), intent(out) :: x1,x2 !output x1, x2
double precision :: pi, s, phi, ra1, ra2 !initialize variables
integer :: i = 0

pi = acos(-1.0) !define pi

do i = 1,d 
	call random_number(ra1)
	call random_number(ra2)
	phi = 2.0*pi*ra1
	s = -delta2*log(1-ra2)
	x1(i) = sqrt(2.0*s)*cos(phi)
	x2(i) = sqrt(2.0*s)*sin(phi)
end do

end subroutine gaussian


subroutine get_parameter
Use ParDataStructure 
use AtomsDataStructure
use NVTDataStructure
!character(len=30), intent(in) :: filename
!
! ----- open input file; read in parameters
	  write(*,*)"Reading parameter file..."
      open( unit=1, file='parameter.par', status='old' )
	  
	  
	  read(1,*) ens
	  read(1,*) config
      read(1,*) N
      read(1,*) maxstep 
	  read(1,*) T_system
      read(1,*) epsk
	  read(1,*) rho
      read(1,*) sigma
      read(1,*) xmass
      read(1,*) timeps
	  
	  if (ens == "NVT") then
		read(1,*)taus
	  end if
!
! ----- factor in energy and force calculation
!
     
! ----- unit eps = K * (J/K) = J
!
      eps  = epsk * boltz
!                                           CONVERSION TO REDUCED UNITS
! ----- reduce well-depth of the LJ potential
!
      epsred = eps / eps
	  T_system = T_system * boltz / eps !convert the temperature of the system to reduced units
	  Tdes = T_system
	  z = 3.0*float(N)*Tdes
	  zeta = 0.0
	  lns = 0.0
	  
	  
!
! ----- reduced mass for a diatomic! molecule in kg
!
	  mu = xmass
      mu   = mu / (avogad * 1000.0 )
!
! ----- reduced time step: 
!   dt = (eps / (mu * sigma^2))^0.5 * (10^-12 s/ps) * (10^+10 Angstrom/m) * timeps
!
      dt   = ((eps / (mu * sigma**2.0))**0.5) * 1.0e-12 * 1.0e+10 * timeps
      dtsq = dt * dt
	  
	  !reduced taus
	  
	  taus   = ((eps / (mu * sigma**2.0))**0.5) * 1.0e-12 * 1.0e+10 * taus
	  
	  
	  time_simulated = dt * maxstep
	  !write(*,*) 'Time simulated (reduced units): ',time_simulated
	  
	  if (ens == "NVT") then
		
		Ms = 3.0*float(N)*T_system*(taus**2.0)
		!Ms = 100000000.0
	  end if
	  
	  !sigma = sigma/eps !reducing sigma?
	  write(*,*)"Got parameters"
      return 
	  

end subroutine



subroutine calc_dist(atom1, atom2, L, dist) !calculate the distance btw two particles vector pointing from atom1 to atom2
use AtomsDataStructure
double precision, dimension(3), intent(out) :: dist !returns the distance as x,y,z vector
integer, intent(in) :: atom1, atom2 !takes in the number of the atoms as argument
double precision, intent(in) :: L !length of the simulation cell, reduced

x1 = coordx(atom1) 
x2 = coordx(atom2)

dist(1) = (x2-x1) - L*anint((x2-x1)/L) !uses pbc


y1 = coordy(atom1) 
y2 = coordy(atom2)

dist(2) = (y2-y1) - L*anint((y2-y1)/L)

z1 = coordz(atom1) 
z2 = coordz(atom2)

dist(3) = (z2-z1) - L*anint((z2-z1)/L)

end subroutine



subroutine nb_list(rn, rc) !!UNUSED
use AtomsDataStructure
use ParDataStructure
double precision, intent(in) :: rn, rc
integer :: i=0, j=0
double precision, dimension(3) :: dist
double precision :: d

do i = 1, N-1
	do j =i+1, N ! only calculating half of the matrix 
		call calc_dist(i,j,L,dist)
		d = (dist(1)**2.0 + dist(2)**2 + dist(3)**2)**(1.0/2.0) !make a scalar from the distance vector
		if (d<rn .AND. d>rc) then
		 ! add to nb list
		
		end if
	end do
end do

end subroutine 



subroutine init_setup_fcc !only this is changed from CC to FCC
Use ParDataStructure
use AtomsDataStructure
double precision, dimension(N) :: temp
double precision, dimension(3) :: dist
double precision :: dcell, delta2, sumvx, sumvy, sumvz, temp2, eta
integer :: ncell, i=0, j=0, k=0, nt=0


rho = rho * avogad/xmass /(10.0**24.0) !rho is read in g/cm3, conversion to 1/A3

ncell = int((N/4.0)**(1.0/3.0)+0.001) !number of subcells
L = (float(N)/rho)**(1.0/3.0) !calculate simulation cell length

Lred = L/sigma !reduce simulation cell

deltar = Lred/float(size(g_hist)) !deltar for g(r)

dcell = Lred / float(ncell) !get subcell length

allocate(coordx(N)) !all needed to be able to put coordinate, velocity, force data
allocate(coordy(N))
allocate(coordz(N))
allocate(coordx_old(N))
allocate(coordy_old(N))
allocate(coordz_old(N))
allocate(velx(N))
allocate(vely(N))
allocate(velz(N))
allocate(Fx(N))
allocate(Fy(N))
allocate(Fz(N))
allocate(forcesx(N,N))
allocate(forcesy(N,N))
allocate(forcesz(N,N))


do i = 1, ncell
	do j = 1, ncell
		do k = 1, ncell 
			nt = nt+1
			coordx(nt) = dcell*float(i-1) !give initial coordinates to the particles
			coordy(nt) = dcell*float(j-1)
			coordz(nt) = dcell*float(k-1)
			
			coordx(nt+1*ncell**3) = dcell*float(i-1)+0.5*dcell !for FCC: each extra atom is shifted half cell side length in 2/3 directions
			coordy(nt+1*ncell**3) = dcell*float(j-1)+0.5*dcell
			coordz(nt+1*ncell**3) = dcell*float(k-1)
			
			coordx(nt+2*ncell**3) = dcell*float(i-1)+0.5*dcell
			coordy(nt+2*ncell**3) = dcell*float(j-1)
			coordz(nt+2*ncell**3) = dcell*float(k-1)+0.5*dcell
			
			coordx(nt+3*ncell**3) = dcell*float(i-1)
			coordy(nt+3*ncell**3) = dcell*float(j-1)+0.5*dcell
			coordz(nt+3*ncell**3) = dcell*float(k-1)+0.5*dcell
			
			
		end do
	end do
end do



delta2 = T_system * mu/mu !calculate deltasquared for gaussian distribution

call gaussian(delta2, N, velx, vely) !give the particles velocities according to a Gaussian dist
call gaussian(delta2, N, velz, temp)


do i =1, N ! shifting the velocities so the COM has 0 velocity
	sumvx = sumvx + velx(i)
	sumvy = sumvy + vely(i)
	sumvz = sumvz + velz(i) 
	
end do


sumvx = sumvx/float(N)
sumvy = sumvy/float(N)
sumvz = sumvz/float(N)

do i =1, N
	velx(i) = velx(i) - sumvx
	vely(i) = vely(i) - sumvy
	velz(i) = velz(i) - sumvz
end do



call calc_temp(temp2)

eta = sqrt(T_system/temp2) !correcting the temperature by scaling

do i =1, N
	velx(i) = velx(i)*eta
	vely(i) = vely(i)*eta
	velz(i) = velz(i)*eta
end do

end subroutine

subroutine init_setup_cc
Use ParDataStructure
use AtomsDataStructure
double precision, dimension(N) :: temp
double precision :: dcell, delta2, sumvx, sumvy, sumvz, temp2, eta
integer :: ncell, i=0, j=0, k=0, nt=0


rho = rho * avogad/xmass /(10.0**24.0) !rho is read in g/cm3, conversion to 1/A3
write(*,*)rho*sigma**3,"density"

L = (N/rho)**(1.0/3.0) !calculate simulation cell length

ncell = int(N**(1.0/3.0)+0.001) !number of subcells
Lred = L/sigma !reduce simulation cell

deltar = Lred/20.0 !deltar for g(r)

dcell = Lred / float(ncell) !get subcell length


allocate(coordx(N)) !all needed to be able to put coordinate, velocity, force data
allocate(coordy(N))
allocate(coordz(N))
allocate(coordx_old(N))
allocate(coordy_old(N))
allocate(coordz_old(N))
allocate(velx(N))
allocate(vely(N))
allocate(velz(N))
allocate(Fx(N))
allocate(Fy(N))
allocate(Fz(N))
!allocate(accx(N))
!allocate(accy(N))
!allocate(accz(N))


do i = 1, ncell
	do j = 1, ncell
		do k = 1, ncell 
			nt = nt+1
			coordx(nt) = dcell*float(i-1) !give initial coordinates to the particles
			coordy(nt) = dcell*float(j-1)
			coordz(nt) = dcell*float(k-1)
		end do
	end do
end do


delta2 = T_system * mu/mu !calculate deltasquared for gaussian distribution

call gaussian(delta2, N, velx, vely) !give the particles velocities according to a Gaussian dist
call gaussian(delta2, N, velz, temp)


do i =1, N ! shifting the velocities so the COM has 0 velocity
	sumvx = sumvx + velx(i)
	sumvy = sumvy + vely(i)
	sumvz = sumvz + velz(i) 
	
end do


sumvx = sumvx/float(N)
sumvy = sumvy/float(N)
sumvz = sumvz/float(N)

do i =1, N
	velx(i) = velx(i) - sumvx
	vely(i) = vely(i) - sumvy
	velz(i) = velz(i) - sumvz
end do


call calc_temp(temp2)

eta = sqrt(T_system/temp2) !correcting the temperature by scaling


do i =1, N
	velx(i) = velx(i)*eta
	vely(i) = vely(i)*eta
	velz(i) = velz(i)*eta
end do

end subroutine




subroutine calc_temp(T) !calculate the temperature in the system from the total kinetic energy
use AtomsDataStructure
use ParDataStructure
double precision, intent(out) :: T
integer :: i=0
double precision :: sumv=0.0

sumv = 0.0
do i =1, N
	sumv = sumv + velx(i)**2.0 + vely(i)**2.0 + velz(i)**2.0 !get the total kinetic energy
end do

T = sumv / 3.0 / float(N) !calculate the temperature, yes N needs to be cast to a double

end subroutine

subroutine calc_pressure(press) 
use AtomsDataStructure
use ParDataStructure
double precision, intent(out) :: press
integer :: i
double precision :: W, p_tail, cutoff

W = 0.0

do i = 1, N-1
	do j =i+1, N ! only calculating half of the matrix 
		W = W + (coordx(j) - coordx(i)) * forcesx(i,j)
		W = W + (coordy(j) - coordy(i)) * forcesy(i,j)
		W = W + (coordz(j) - coordz(i)) * forcesz(i,j)
	end do
end do

press = W
do i = 1, N

press = press + (velx(i)**2)+ (vely(i)**2) + (velz(i)**2)

end do

press = 1.0/3.0 * press /(Lred**3.0)

if (isTailOn .eqv. .True.) then

cutoff = 1.0/2.5 !this is actually 1/rc here
p_tail = 32.0/9.0 * pi*rho**2.0*(sigma/eps)**3.0*(cutoff**9.0-1.5*cutoff**3.0)
press = press + p_tail

end if

end subroutine




subroutine ener_force_NVE
use AtomsDataStructure
use ParDataStructure
double precision, dimension(3) :: dist, dist_test
double precision :: var1,var2,distance,force,cutoff
integer :: ii
!double precision, dimension(N,N) :: forcesx, forcesy, forcesz

forcesx = 0.0 !zero out the forces, these include Fij and Fji
forcesy = 0.0
forcesz = 0.0


! ----- calculate energy and force

Fx = 0.0 !these the total forces acting on atom i
Fy = 0.0
Fz = 0.0

Epot = 0.0 !reset Epot before counting everything

cutoff = 2.5*sigma

do i = 1, N-1 !the problem is probably somwhere here since everything depending on the force is messed up
	do j =i+1, N ! only calculating half of the matrix 
		!calculate forces between pairs
		call calc_dist(i,j,Lred,dist) !get the distance vector btw the two atoms, does use Lred for cell length as it should
		distance = sqrt(dist(1)**2+dist(2)**2+dist(3)**2)
		if (imd >10000) then
		ii = anint(distance/deltar) + 1
		g_hist(ii) = g_hist(ii) + 2 !add to the histogram, maybe move this to analysis
		end if 
		!write(*,*)distance,'distance btw the atoms'
		if (isTailOn .eqv. .True.) then!if I want to implement cutoff then this happens
		
			if (distance < cutoff) then
				var1 = (1.0/distance)**(13.0) !calculate the force from the Lennard-Jones potential
				var2 = -(1.0/distance)**(7.0)
				force = -48.0*(var1+0.5*var2)
		
				var1 = (1.0/distance)**(12.0) !add the potential from the force to the total potential
				var2 = (1.0/distance)**(6.0) 
				Epot = Epot + 4.0*(var1-var2)
		

				forcesx(i,j) = dist(1)/distance * force
				forcesy(i,j) = dist(2)/distance * force
				forcesz(i,j) = dist(3)/distance * force

				forcesx(j,i) = -forcesx(i,j)  !give the other half of the pair 
				forcesy(j,i) = -forcesy(i,j)
				forcesz(j,i) = -forcesz(i,j) 
		
			end if
		else
		
		var1 = (1.0/distance)**(13.0) !calculate the force from the Lennard-Jones potential
		var2 = -(1.0/distance)**(7.0)
		force = -48.0*(var1+0.5*var2)
		
		var1 = (1.0/distance)**(12.0) !add the potential from the force to the total potential
		var2 = (1.0/distance)**(6.0) 
		Epot = Epot + 4.0*(var1-var2)
		

		forcesx(i,j) = dist(1)/distance * force
		forcesy(i,j) = dist(2)/distance * force
		forcesz(i,j) = dist(3)/distance * force

		forcesx(j,i) = -forcesx(i,j)  !give the other half of the pair 
		forcesy(j,i) = -forcesy(i,j)
		forcesz(j,i) = -forcesz(i,j) 
		
		end if
	end do
end do	

do i = 1,N !sum up the forces acting on the atoms
	do j = 1,N
	Fx(i) = Fx(i) + forcesx(i,j)
	Fy(i) = Fy(i) + forcesy(i,j)
	Fz(i) = Fz(i) + forcesz(i,j)
	end do
end do

end subroutine

subroutine ener_force_NVT !NOT USED
use AtomsDataStructure
use ParDataStructure
use NVTDataStructure
double precision, dimension(3) :: dist, dist_test
double precision :: var1,var2,distance,force
!double precision, dimension(N,N) :: forcesx, forcesy, forcesz


forcesx = 0.0 !zero out the forces, these include Fij and Fji
forcesy = 0.0
forcesz = 0.0

!if (.not. imd == 1) then
!	do i=1, N !output the old acceleration for the velocity Verlet before zeroing out the force
!		accx(i) = Fx(i) - zeta *velx(i)
!		accy(i) = Fy(i) - zeta *vely(i)
!		accz(i) = Fz(i) - zeta *velz(i)
!	end do
!end if

! ----- calculate energy and force

Fx = 0.0 !these the total forces acting on atom i
Fy = 0.0
Fz = 0.0

Epot = 0.0 !reset Epot before counting everything

do i = 1, N-1 
	do j =i+1, N ! only calculating half of the matrix 
		!calculate forces between pairs
		call calc_dist(i,j,Lred,dist) !get the distance vector btw the two atoms, does use Lred for cell length as it should
		distance = sqrt(dist(1)**2+dist(2)**2+dist(3)**2)
		!write(*,*)distance,'distance btw the atoms'
		
		var1 = (1.0/distance)**(13.0) !calculate the force from the Lennard-Jones potential
		var2 = -(1.0/distance)**(7.0)
		force = -48.0*(var1+0.5*var2)
		
		var1 = (1.0/distance)**(12.0) !add the potential from the force to the total potential
		var2 = (1.0/distance)**(6.0) 
		Epot = Epot + 4.0*(var1-var2)
		

		forcesx(i,j) = dist(1)/distance * force
		forcesy(i,j) = dist(2)/distance * force
		forcesz(i,j) = dist(3)/distance * force

		forcesx(j,i) = -forcesx(i,j)  !give the other half of the pair 
		forcesy(j,i) = -forcesy(i,j)
		forcesz(j,i) = -forcesz(i,j) 
		
		if (imd > 10000) then
		!g(r)
		
		end if
	end do
end do	

do i = 1,N !sum up the forces acting on the atoms
	do j = 1,N
	Fx(i) = Fx(i) + forcesx(i,j)
	Fy(i) = Fy(i) + forcesy(i,j)
	Fz(i) = Fz(i) + forcesz(i,j)
	end do
end do

!if (imd == 1 ) then !just for the first step 
!	do i=1, N 
!		accx(i) = Fx(i) 
!		accy(i) = Fy(i) 
!		accz(i) = Fz(i) 
!	end do
!end if

end subroutine




subroutine update_NVE
use AtomsDataStructure
use ParDataStructure
double precision :: temp
double precision, dimension(3) :: old_old

! ----- update velocity and coordinate
!

Ekin = 0.0 !zero out the total kinetic energy

do i = 1,N
      if( imd == 1 ) then !in the first step there needs to be a -1st coordinate calculated
	  
	  !VERLET
		old_old(1)   = coordx(i) - dt * velx(i) + 0.5 * Fx(i) * dtsq
		old_old(2)   = coordy(i) - dt * vely(i) + 0.5 * Fy(i) * dtsq
		old_old(3)   = coordz(i) - dt * velz(i) + 0.5 * Fz(i) * dtsq
		
		temp = 2.0 * coordx(i) - old_old(1) + Fx(i)*dtsq !calc the new coordinate
		coordx_old(i) = coordx(i) !put the old coordinate to the right variable
		coordx(i) = temp !put the new coordinate to the right variable
		
		temp = 2.0 * coordy(i) - old_old(2) + Fy(i)*dtsq
		coordy_old(i) = coordy(i)
		coordy(i) = temp
		
		temp = 2.0 * coordz(i) - old_old(3) + Fz(i)*dtsq
		coordz_old(i) = coordz(i)
		coordz(i) = temp

      else
		! VERLET 
		
		temp = 2.0 * coordx(i) - coordx_old(i) + Fx(i)*dtsq	!calculate the new coordinate
		old_old(1) = coordx_old(i) !put the old coord to the very old coord
		coordx_old(i) = coordx(i) !put the old coord to the right variable
		coordx(i) = temp !put the new coord to the right variable
		
		temp = 2.0 * coordy(i) - coordy_old(i) + Fy(i)*dtsq		
		old_old(2) = coordy_old(i) 
		coordy_old(i) = coordy(i)
		coordy(i) = temp
		
		temp = 2.0 * coordz(i) - coordz_old(i) + Fz(i)*dtsq		
		old_old(3) = coordz_old(i) 
		coordz_old(i) = coordz(i)
		coordz(i) = temp
		
      end if

	velx(i) = 0.5/dt * (coordx(i)-old_old(1)) !calculate the velocity
	vely(i) = 0.5/dt * (coordy(i)-old_old(2))
	velz(i) = 0.5/dt * (coordz(i)-old_old(3))
		
	Ekin = Ekin + 0.5 * (velx(i)**2 + vely(i)**2 + velz(i)**2) !sum up the kinetic energy in the system)

end do

hamilt = Epot + Ekin !this is okay bc the kinetic energy is calculated after checking the potential before

call calc_temp(T_system)
!call calc_pressure(p_system)
write(10,'(6f12.6)')(imd-1.0)*dt,Epot,Ekin,hamilt,T_system
h_t(imd) = hamilt
end subroutine

subroutine update_NVT
use AtomsDataStructure
use ParDataStructure
use NVTDataStructure
use accelerations

double precision :: temp, sumv, sumva, zetadot, denom
double precision, dimension(3) :: old_old
integer :: i

! ----- update velocity and coordinate
!

Ekin = 0.0 !zero out the total kinetic energy

!if (.not. imd == 1) then
accx=Fx - zeta *velx
accy=Fy - zeta *vely
accz=Fz - zeta *velz

!end if

!write(*,*)imd,T_system,"T"
call calc_temp(T_system)
sumv = 3.0*float(N)*T_system !T_system at t
!write(*,*)imd,sumv,"sumv"
!write(*,*)imd,z,"z"
sumva = 0.0
do i = 1,N
	!sumv = sumv + velx(i)**2.0 + vely(i)**2.0 + velz(i)**2.0
	sumva = sumva + velx(i)*accx(i) + vely(i)*accy(i) + velz(i)*accz(i) !sumva at t
end do

zetadot = (1.0/Ms) * (sumv - z) !calculate zetadot at t
!update lns for the Hamiltonian
lns = lns + zeta*dt + zetadot*dtsq/2.0 !update ln s from t to t+dt


!write(*,*)imd,zetadot,"zetadot"

!update the friction coefficient
zeta = zeta + (dt*zetadot) + (dtsq/Ms * sumva) !update zeta from t to t+dt


!write(*,*)imd,zeta,"zeta"

denom = 1.0 + zeta*dt/2.0

!write(*,*)imd,denom



do i = 1,N

! for some godforsaken reason the algorithm doesn't work if I use the velocity verlet to move the coordinates. the regular Verlet is also O(n^4) error so I will keep using rhe old verlet.
      !if( imd == 1 ) then !in the first step there needs to be a -1st coordinate calculated
	  
	  !VERLET 
		!old_old(1)   = coordx(i) - dt * velx(i) + 0.5 * Fx(i) * dtsq
		!old_old(2)   = coordy(i) - dt * vely(i) + 0.5 * Fy(i) * dtsq
		!old_old(3)   = coordz(i) - dt * velz(i) + 0.5 * Fz(i) * dtsq
		
		!temp = 2.0 * coordx(i) - old_old(1) + Fx(i)*dtsq !calc the new coordinate
		!coordx_old(i) = coordx(i) !put the old coordinate to the right variable
		!coordx(i) = temp !put the new coordinate to the right variable
		
		!temp = 2.0 * coordy(i) - old_old(2) + Fy(i)*dtsq
		!coordy_old(i) = coordy(i)
		!coordy(i) = temp
		
		!temp = 2.0 * coordz(i) - old_old(3) + Fz(i)*dtsq
		!coordz_old(i) = coordz(i)
		!coordz(i) = temp

     ! else
		! VERLET 
		
		!temp = 2.0 * coordx(i) - coordx_old(i) + Fx(i)*dtsq	!calculate the new coordinate
		!old_old(1) = coordx_old(i) !put the old coord to the very old coord
		!coordx_old(i) = coordx(i) !put the old coord to the right variable
		!coordx(i) = temp !put the new coord to the right variable
		
		!temp = 2.0 * coordy(i) - coordy_old(i) + Fy(i)*dtsq		
		!old_old(2) = coordy_old(i) 
		!coordy_old(i) = coordy(i)
		!coordy(i) = temp
		
	!	temp = 2.0 * coordz(i) - coordz_old(i) + Fz(i)*dtsq		
	!	old_old(3) = coordz_old(i) 
	!	coordz_old(i) = coordz(i)
	!	coordz(i) = temp
		
     ! end if

	coordx(i) = coordx(i) + velx(i)*dt +accx(i)*(dt**2.0)/2.0	!calculate the new coordinate
		
	coordy(i) = coordy(i) + vely(i)*dt +accy(i)*(dt**2.0)/2.0	!calculate the new coordinate 
		
	coordz(i) = coordz(i) + velz(i)*dt +accz(i)*(dt**2.0)/2.0	!calculate the new coordinate

		

	
	
	
	

end do


call ener_force_NVE !move F from t do t+dt

do i = 1,N

	temp = velx(i) + (dt*accx(i))/2.0 + (Fx(i)*dt)/2.0 !velocity at t+dt
	velx(i) = temp/denom
	temp = vely(i) + (dt*accy(i))/2.0 + (Fy(i)*dt)/2.0
	vely(i) = temp/denom
	temp = velz(i) + (dt*accz(i))/2.0 + (Fz(i)*dt)/2.0
	velz(i) = temp/denom
		
	Ekin = Ekin + 0.5 * (velx(i)**2 + vely(i)**2 + velz(i)**2) !sum up the kinetic energy in the system)

end do

!write(*,*)imd,"Epot",Epot,"Ekin",Ekin
!write(*,*)imd,"lns",lns,"zeta",zeta

hamilt = Epot + Ekin + z*lns +Ms*(zeta**2.0)/2.0
!hamilt = Epot + Ekin
!this is okay bc the kinetic energy is calculated after checking the potential before

call calc_temp(T_system)
call calc_pressure(p_system)
write(10,'(6f12.6)')(imd)*dt,Epot,Ekin,hamilt,T_system,zeta
h_t(imd) = hamilt
comp = p_system*(Lred)**3.0/(N*T_system) !compression factor
z_t(imd) = comp
end subroutine



subroutine analysis(slope)
use AtomsDataStructure
use ParDataStructure
use NVTDataStructure
integer :: i, ii
double precision, dimension(imd) :: h_t_red, times 
double precision, intent(out) :: slope
double precision :: intercept

h_t_red = 0.0

do i = 1,imd
	h_t_red(i) = h_t(i) !cuts off the 0s from the array containing the Hamiltonian history
	times(i) = (i-1)*dt
end do

call linreg(times,h_t_red,slope)



end subroutine

subroutine linreg(Lx,Ly,m)
use ParDataStructure
double precision, intent(in), dimension(imd)  :: Lx, Ly
double precision, intent(out) :: m
double precision :: sum_x, sum_y, sum_xx, sum_xy


sum_x  = sum(Lx)
sum_y  = sum(Ly)

sum_xx = sum(Lx * Lx)
sum_xy = sum(Lx * Ly)
m = (sum_xy * imd - sum_y * sum_x) &
	& / (imd * sum_xx - sum_x * sum_x)

end subroutine

subroutine goodbye
double precision:: text_index

call random_number(text_index)
write(*,*)" "

if (text_index<0.1) then
	write(*,*)'MOOMIN reminds you: "FORTRAN is not a flower but a weed - it is hardy,&
	&occasionally blooms, and grows in every computer" - Alan J. Perlis'
else if (text_index<0.2) then
	write(*,*)'MOOMIN reminds you: "Recall the first American space probe to Venus, reportedly lost &
	&because Fortran cannot recognize a missing comma in a DO statement… " - C. A. R. Hoare'
else if (text_index<0.3) then
	write(*,*)'MOOMIN reminds you: "The determined Real Programmer &
	&can write FORTRAN programs in any language." - Ed Post'
else if (text_index<0.4) then
	write(*,*)'MOOMIN reminds you: "FORTRAN was the language of choice for the &
	&same reason that three-legged races are popular." - Ken Thompson'
else if (text_index<0.5) then
	write(*,*)'MOOMIN reminds you: "God is real, unless declared otherwise." - the depths of Stackoverflow'
else if (text_index<0.6) then
	write(*,*)'MOOMIN reminds you: "We are nothing but the sand that fills the cracks, &
	&a soul in denial, a body ashamed" - Moonspell'
else if (text_index<0.7) then
	write(*,*)'MOOMIN reminds you: "So many stones out there to turn, &
	&and yet we choose to ignore" - Moonspell '
else if (text_index<0.8) then
	write(*,*)'MOOMIN reminds you: "Is it right to indulge on an ecstasy &
	&of creating a god that sees what I see?" - Moonspell'
else if (text_index<0.9) then
	write(*,*)'MOOMIN reminds you: "Slow down. God can hear you." - Moonspell'
else
	write(*,*)'MOOMIN reminds you: "STOP programming computers. Sand was never &
	&meant to think. This is very cruel to rocks." - r/SurrealMemes'
end if

end subroutine


program moomin
use AtomsDataStructure
use ParDataStructure
use NVTDataStructure
use accelerations
double precision :: text_index, drift, temp
double precision, allocatable, dimension(:) :: g_hist_out
integer :: i, ii, imd_i, j

pi = acos(-1.0) !define pi


call random_number(text_index)

write(*,*)":-) MOOMIN - ./moomin.out, 2022-Ubuntu (-:"

if (text_index<0.2) then
	write(*,*)"MOOMIN stands for: MOnomOlecular Material INtegrator"
else if (text_index<0.4) then
	write(*,*)"MOOMIN stands for: runtime Magnitude: Order Of MIllioN seconds"
else if (text_index<0.6) then
	write(*,*)"MOOMIN stands for: Making Obvious acrOnyMs Is Not easy"
else if (text_index<0.8) then
	write(*,*)"MOOMIN stands for: MOOnspell Motivated INnnovation"
else
	write(*,*)"MOOMIN stands for: MOlekylær dynamik simulering Og MachINe learning"
end if

write(*,*)" "
write(*,*)"MOOMIN is written by:"

write(*,*)"Alíz Hanga Lelik"

write(*,*)"with the help of Günther H.J. Peters"

write(*,*)"for the course 26255 Molecular dynamics and machine learning at DTU" 
write(*,*)" "

write(*,*)"MOOMIN is free software, unless you want to pay for it.&
& But why would you? Go use GROMACS for any serious science."
write(*,*)" "


!open output file
open( unit=10, file='output.dat', status='unknown' )
open( unit=11, file='output_gr.dat', status='unknown' )

call get_parameter


if (config == "FCC") then
call init_setup_fcc
else
call init_setup_cc
end if
write(*,*)"Calculating..."

allocate(accx(N))
allocate(accy(N))
allocate(accz(N))


allocate(h_t(maxstep))
allocate(z_t(maxstep))
h_t = 0.0
z_t = 0.0

g_hist = 0
allocate(g_hist_out(size(g_hist)))
g_hist = 0.0

isTailOn = .False.


temp = float(maxstep)/100.0
ii = int(temp-1) !divide the main do loop into two so the analysis is only performed every 100 steps

if (ens == "NVE") then
	do imd_i=1, ii
		do i = 1, 100
		imd = (imd_i-1)*100+i

		call ener_force_NVE

		call update_NVE
		end do
		
		call analysis(drift)
		if ((drift>0.1) .or. (drift<-0.1)) then
		write(*,*)"Energy drift out of bounds. Stopping program."
		call goodbye
		stop
		end if
		
		WRITE(*,101, ADVANCE='NO')achar(13),imd, maxstep
		101     FORMAT( a , ' Step number : ',i7,' out of ',i7)
	end do

	ii = maxstep-imd

	do i=1,ii
		imd = imd +1
		call ener_force_NVE
		call update_NVE
	end do
	
elseif (ens == "NVT") then
	imd = 1
	call ener_force_NVE
	!write(*,*)"ener_force"
	do imd_i=1, ii
		do i = 1, 100
		imd = (imd_i-1)*100+i
		
		
		!if (all(Fx .eq. accx)) then
		!	write(*,*)imd,"accelerations correctly casted"
		!end if
		!call ener_force_NVE
		
		!if (imd == 1 ) then !just for the first step 
			
		!	accx = Fx
		!	accy = Fy
		!	accz = Fz
		!end if
		
		call update_NVT
		end do
		
		call analysis(drift)
		z_t(imd_i) = comp
		
		WRITE(*,101, ADVANCE='NO')achar(13),imd, maxstep
	end do

	ii = maxstep-imd

	do i=1,ii
		imd = imd +1
		if (.not. imd == 1) then
			accx=Fx - zeta *velx(i)
			accy=Fy - zeta *vely(i)
			accz=Fz - zeta *velz(i)
		end if
		
		call ener_force_NVE
		call update_NVT
	end do
else
!other ensembles
end if
WRITE(*,101, ADVANCE='NO')achar(13),imd, maxstep


write(*,*)" "
write(*,*)" "
write(*,*)"Total energy of the system: ",hamilt
write(*,*)"Total energy drift: ",drift
write(*,*)"Final system temperature: ",T_system*epsk," K"
comp = sum(z_t)/max(1,size(z_t))
write(*,*)"Average compression factor: ",comp
write(*,*)" "

!write(*,*)g_hist
do i=1,size(g_hist)

	temp = float(N)*float(maxstep-10000)*rho*(4.0*pi)*(deltar**3.0)*(((i+1)**3.0-i**3.0)/3.0)
	!write(*,*)temp
	g_hist_out(i) = float(g_hist(i))/temp
	write(11,*)float((i-1))*deltar,g_hist_out(i)
end do


call goodbye

end program