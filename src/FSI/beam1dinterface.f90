! example 6.4 on p. 237 of Reddy
subroutine beam1d(LociRank, LociNNodes, LociDofFree, CSDnodes, plungingType, plungeAmp, flapAmp, freq, chord, rho, thickness, Young, CSDx, CSDxdot, CSDxddot, CSDfpre, CSDForce, dt, sTime, timestep, xNew, xdotNew, xddotNew, FpreNew, disp) 
use mesh
use stiff
use gauss
!use eigen
!use matrix_inverter
implicit none

! declarations, passed from CFD
integer :: LociRank ! number of nodes
integer :: LociNNodes ! number of nodes
integer :: LociDofFree ! = dofFreetotal degree of freedom 
integer :: plungingType ! 0 = switch off, 1 = sin, 2 = cos, 3 = 1-cos
real, dimension(LociNNodes) :: CSDnodes
integer, dimension(:,:), allocatable :: bcIndex ! indices of the imposed boundary conditions: bcIndex(dofNode, node)
real :: plungeAmp, flapAmp
real  :: freq ! f, not omega
real  :: chord
real  :: rho
real  :: thickness
real  :: Young
real  :: dt ! time step size
real  :: sTime ! simulation time
integer  :: timestep ! current time step #
real, dimension(LociDofFree) :: CSDx
real, dimension(LociDofFree) :: CSDxdot
real, dimension(LociDofFree) :: CSDxddot
real, dimension(LociDofFree) :: CSDFpre
real, dimension(LociNNodes) :: CSDForce

! declarations, passed to CFD
real, dimension(LociDofFree) :: xNew
real, dimension(LociDofFree) :: xdotNew
real, dimension(LociDofFree) :: xddotNew
real, dimension(LociDofFree) :: FpreNew
real, dimension(LociNNodes) :: disp

! declarations, inside beam1d
real :: area, inertia
!integer :: nElements ! number of elements
!integer :: dofElement ! degree of freedom per element
!integer :: nodePerElement ! number of nodes per element
!integer :: nBC ! number of Boundary conditions -> 2 at the root
!integer :: dofNode ! number of dof per node
!integer :: dofFree ! number of free dof free -> dofNode*nNodes - nBC
!integer :: LociDofFree ! =dofFree from Loci-STREAM: number of free dof free -> dofNode*nNodes - nBC
integer :: n, e, i, t ! indices nodes/elements
integer :: nTimeSteps ! number of timesteps, here always 1
integer :: ptrOut=84, ptrOutTip=85
integer :: outFreq
!integer :: rank

real,dimension(LociDofFree) :: x
real,dimension(LociDofFree) :: xdot
real,dimension(LociDofFree) :: xddot
real,dimension(LociDofFree) :: Fpre

real :: PI
!real :: alpha, gammaNM ! Newmark parameters: gammaNM -> gamma from NewMark  
real :: alpha,beta, gammaT ! Newmark parameters: gammaNM -> gamma from NewMark  
real :: q ! distribution
!real :: a1, a2, a3, a4, a5 ! Newmark parameters
real :: c1, c2, c3, alphaP1 ! Newmark parameters
real, dimension(:,:), allocatable :: K ! global stiffness matrix
real, dimension(:,:), allocatable :: M ! global mass matrix
real, dimension(:), allocatable :: Pn, F ! element distributed and point force vector and global force vector
real, dimension(:,:), allocatable :: A, v !  M^-1 F , eigenvectors
real, dimension(:), allocatable :: d ! eigenvalues
!integer :: nrot ! number of Jacobi rotations needed to compute eigensystem

!real, parameter :: E = ! Young's modulus 
!real, parameter :: I = ! First moment of inertia
real :: EIdefault, rhoAdefault 

real, dimension(:,:), allocatable :: sol ! solution array: 1-> dof, 2-> nodes

character(len=80) :: outputFile, outputFileTip, temp

xNew=0.0d0
xdotNew=0.0d0
xddotNew=0.0d0
FpreNew=0.0d0

rank = LociRank 
nNodes = LociNNodes
dofNode = 2 
nBC = 2
dofFree = nNodes*dofNode - nBC
nElements = nNodes - 1
nodePerElement = 2
dofElement = dofNode * nodePerElement

PI = 2. * acos(0.0d0)
! initialization


nTimeSteps = 1
outFreq = 1

!sTime = 0.0d0
!alpha = 0.5d0
!gammaNM = 0.5d0
alpha = -0.1; beta = 0.3025; gammaT = 0.6 ! see Hilber1

! derived

q = 0.0d0
area = thickness * chord
inertia = chord * thickness**3 / 12.0d0
EIdefault = Young * inertia
rhoAdefault = rho * area

if (rank==0) print*, 'nNodes = ', nNodes
if (rank==0) print*, 'timestep = ', timestep
if (rank==0) print*, 'dt = ', dt
if (rank==0) print*, 'stime = ', stime
if (rank==0) print*, 'EI = ', EIdefault
if (rank==0) print*, 'rhoA = ', rhoAdefault
if (rank==0) print*, 'plungeAmp = ', plungeAmp
if (rank==0) print*, 'flapAmp = ', flapAmp
if (rank==0) print*, 'freq = ', freq
if (rank==0) then 
      if (dofFree == LociDofFree) print*, 'dofFree correct ', dofFree
endif

allocate( nodes(nNodes) )
nodes=CSDnodes
x = CSDx
xdot = CSDxdot
xddot = CSDxddot
Fpre = CSDFpre

!if (rank==0) print*, '[I] nodes', nodes
!if (rank==0) print*, '[I] x', x
!if (rank==0) print*, '[I] xdot', xdot
!if (rank==0) print*, '[I] xddot', xddot



outputFileTip = 'tipoutput.dat'



! allocate
allocate( sol(dofNode, nNodes), K(dofFree, dofFree), Pn(nNodes), F(dofFree), M(dofFree, dofFree), bcIndex(2,nBC))
K=0.0d0; Pn=0.0d0 ; F=0.0d0; M=0.0d0; sol=0.0d0; bcIndex=0

! allocate related to modal analysis
!allocate( A(dofFree, dofFree), v(dofFree, dofFree), d(dofFree))

! boundary conditions (/ dof, node /)
! bcIndex(:,1) = (/1, 1/); sol(1, 1) = plungeAmp
! bcIndex(:,2) = (/2, 1/); sol(2, 1) = flapAmp
bcIndex(:,1) = (/1, 1/); sol(1, 1) = 0.0d0
bcIndex(:,2) = (/2, 1/); sol(2, 1) = 0.0d0
!bcIndex(:,3) = (/2, nNodes/); sol(2, nNodes) = 0.0d0
if (rank==0) print*, '[I] EulerBeam1D: Boundary conditions assigned'

! generate mesh
call assemble_mesh(EIdefault, rhoAdefault)
if (rank==0) print*, '[I] EulerBeam1D: Mesh succesfully assembled'

if (timeStep == 1) then
	! time integration
	! initial conditions
!	x = 0.0d0 ; xdot=0.0d0 
	do n=1,nNodes
	    if (nodeToEquation(1,n) /= -1) then  
				x(nodeToEquation(1, n)) = 0.0d0
!				x(nodeToEquation(1, n)) = plungeAmp
	    end if
	    if (nodeToEquation(2,n) /= -1) then  
				x(nodeToEquation(2, n)) = 0.0d0
!				x(nodeToEquation(2, n)) = flapAmp
	    end if 
	end do
	!x(nodeToEquation(1,nNodes)) = 0.2146
	!xdot = 0.0d0
	if (rank==0) print*, '[I] EulerBeam1D: Initial conditions applied 0'
	call assemble_matrices(sol, q, Pn, M, K, F)
	if (rank==0) print*, '[I] EulerBeam1D: Initial conditions applied 1'


	xddot = F - matmul(K,x)
	call lin_solve(M, xddot)

	! Point load assumed to be zero for the approximation of the accelerations
	!Pn(nNodes) = 60.0e3

	if (rank==0) print*, '[I] EulerBeam1D: Initial conditions applied 2'
	    ! output the initial conditions
	do n=1,nNodes
		do i=1,dofNode
		    if ( nodeToEquation(i, n) /= -1 ) then
			    sol(i,n) = x( nodeToEquation(i,n) )
		    end if
		end do
	end do
	if (rank==0) then
	    write(temp, '(I5.5)') 0
	    outputFile = 'output'//trim(temp)//'.dat'
	    print*, 'Writing the field to '//outputFile
	    open(unit=ptrOut, file=trim(outputFile), action='write')

	    write(ptrOut, '(A)') 'variables = "x", "w", "-theta"'
	    write(ptrOut, '(A, I5, A)') 'zone i=', nNodes, ' datapacking=block'

	    write(ptrOut, *) (nodes(n), n=1,nNodes)
	    write(ptrOut, *) (sol(1, n), n=1,nNodes)
	    write(ptrOut, *) (sol(2, n), n=1,nNodes)
	    close(ptrOut)
	!    print*, sTime, x(nodeToEquation(1, elements(2, nElements)))
	    
	    print*, 'Writing the field to '//outputFileTip
	    open(unit=ptrOutTip, file=trim(outputFileTip), action='write')
	    write(ptrOutTip, '(I5, 2E15.7)') 0, sol(1, 1), sol(1, nNodes)
	    close(ptrOutTip)
	endif
endif

Pn=CSDForce ! In the first time step, for the approximation of the acceleration, the point force is assumed to be zero

if (rank==0) print*, '[I] EulerBeam1D: Time integration starts'
!do t=timestep,(timestep+nTimeSteps-1)
t = timestep

!    a1 = alpha * dt; a2 = (1.0d0-alpha) * dt; a3 = 2.0d0 / (gammaNM * dt**2); a4 = 2.0d0 / (gammaNM * dt); a5 = 1.0d0/gammaNM - 1.0d0

    ! boundary conditions
    if (plungingType==1) then
      sol(1,1) = plungeAmp * ( sin(2 * PI * freq * sTime) )
      if (rank==0) print*, '[I] EulerBeam1D: Sin motion'
!    elseif (flapptingType==2) then
!      sol(1,1) = plungeAmp * ( cos(2 * PI * freq * sTime) )
    elseif (plungingType==3) then
      sol(1,1) = plungeAmp * ( -1.0d0 + cos(2 * PI * freq * sTime) )
      if (rank==0) print*, '[I] EulerBeam1D: Cos - 1 motion'
    endif

     sol(2,1) = flapAmp * ( sin(2 * PI * freq * sTime) )
!     sol(1,1) = plungeAmp * ( -1.0d0 + cos(2 * PI * freq * sTime) )
!     sol(2,1) = flapAmp * ( -1.0d0 + cos(2 * PI * freq * sTime) )
!    sol(nNodes,1) = 0.0d0
   
    call assemble_matrices(sol, q, Pn, M, K, F)

! update force to CFD
fpreNew = F 
    
!     K = K + a3 * M
!     xNew = a3 * x + a4 * xdot + a5 * xddot
!     xNew = matmul(M, xNew)
!     xNew = F + xNew

    alphaP1 = alpha + 1.d0
    F = F*alphaP1 - Fpre * alpha
    c1 = dt**2 * beta
    c2 = dt**2 * (0.5d0 - beta)
    c3 = 1.d0 - gammaT
    M = M + alphaP1 * c1 * K
    xddotNew = alphaP1 * ( x + dt * xdot + c2 * xddot) 
    xddotNew = matmul(K, xddotNew)
    xddotNew = F - xddotNew + alpha * matmul(K, x)

    ! solve
    if (rank==0) print*, '[I] EulerBeam1D: Solving the equations'
    !if (rank==0) print*, '[I] EulerBeam1D: A = ', K
!    call lin_solve(K, xNew)
  
    call lin_solve(M, xddotNew)
    
    if (rank==0) print*, '[I] EulerBeam1D: Equations solved'	

    ! velocity and accelerations for the new time step
!     xddotNew = a3 * ( xNew - x) - a4 * xdot - a5 * xddot
!     xdotNew = xdot + a2 * xddot + a1 * xddotNew
    xNew = x + dt * xdot + c2 * xddot + c1 * xddotNew
    xdotNew = xdot + dt * ( c3 * xddot + gammaT * xddotNew)

		if (rank==0) print*, '[I] EulerBeam1D: New velocities and accelerations computed'
		
    ! update
    !xddot = xddotNew
    !xdot = xdotNew
    !x = xNew
    	
    do n=1,nNodes
        do i=1,dofNode
            if ( nodeToEquation(i, n) /= -1 ) then
                    sol(i,n) = xNew( nodeToEquation(i,n) )
            end if
        end do        
    end do
		
		disp = sol(1, :)
		
		!print*, '[I] EulerBeam1D rank = ', rank, 'disp = ', disp

		if (rank==0) print*, '[I] EulerBeam1D: Displacements to CFD computed', rank, '=rank'
		!print*, '[I] EulerBeam1D: Displacements to CFD computed', rank, '=rank'
			
    ! output quarter and 3/4 span displacements

!        write(ptrOut, '(A)') 'variables = "x", "w", "-theta"'
!        write(ptrOut, '(A, I5, A)') 'zone i=', nNodes, ' datapacking=block'
		if (rank==0) then
        open(unit=ptrOutTip, file=trim(outputFileTip), action='write', access='APPEND')
        write(ptrOutTip, '(I5, 2E15.7)') t, sol(1, 1),  sol(1, nNodes)
        close(ptrOutTip)
        print*, '[I] EulerBeam1D: Writing the tip outputs'

    ! output node positions
	    if (mod(t, outFreq) == 0) then        
	        write(temp, '(I5.5)') t
	        outputFile = 'output'//trim(temp)//'.dat'
	        print*, 'Writing the field to '//outputFile
	        open(unit=ptrOut, file=trim(outputFile), action='write')
	
	        write(ptrOut, '(A)') 'variables = "x", "w", "-theta"'
	        print*, 'nNodes = ', nNodes
	        write(ptrOut, '(A, I5, A)') 'zone i=', nNodes, ' datapacking=block'
	
	        write(ptrOut, *) (nodes(n), n=1,nNodes)
	        write(ptrOut, *) (sol(1, n), n=1,nNodes)
	        write(ptrOut, *) (sol(2, n), n=1,nNodes)
	        close(ptrOut)
	    end if
    
		end if
    
!    print*, sTime, nodeToEquation(1,elements(2,nElements)), x(nodeToEquation(1, elements(2, nElements)))
    
!end do ! time step

! update

if (rank==0) print*, '[I] EulerBeam1D: Completed'

call deallocate_mesh

deallocate(sol, K, F, Pn, M,  bcIndex)
       
end subroutine beam1d

