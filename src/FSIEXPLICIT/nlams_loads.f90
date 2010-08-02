!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_loads. Satish Chimakurthi. 22Aug2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   Read forces and moments applied at nodal location from an input fule 
!
!-> Reference:
!
!->Subroutines:
!
!     |-nlams_nodal_loads
!
!->Remarks:
! 1) Reads forces (3) and moments (3) applied at nodal locations from an input file
! called loads.dat
! 2) Reads nodal boundary conditions from a file called bc.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_loads
! Global variables
 integer,save                   ::NumForceNodes     !Number of nodes which have force/moment BC
 integer,save                   ::NumBC             !Number of boundary conditions
 integer,save,allocatable       ::bcdof(:)
 real(8),save,allocatable       ::bcval(:)
 real(8),save,allocatable       ::RHS(:,:)
!end module nlams_loads

contains
subroutine nlams_nodal_loads
 use mesh_connect_vars,only:NumNodes
! use nlams_loads
 use LociSTREAMData, only: LociSTREAMNumBC,LociSTREAMBCdof,LociSTREAMBCval

 implicit none
 integer:: i

!    allocate(RHS(NumNodes*6,1)); RHS = 0.d0 ;
!    NumForceNodes = 0
!    NumBc = LociSTREAMNumBC;
!    allocate(bcdof(NumBC)); bcdof = 0 ;
!    allocate(bcval(NumBC)); bcval = 0.d0 ;
!    do i=1,NumBc
!      bcdof(i) = LociSTREAMBCdof(i)
!      bcval(i) = LociSTREAMBCval(i)
!    enddo


!Local variables
! integer:: io,np,num,i,j
! real(8),save,allocatable::NodalLoads(:,:)
! integer,save,allocatable::NodeNumb(:)
! Switch off this function since no applied force on the wing
! Input element coordinates
!  open(1022,file='loads.dat',action='read')
!  np=0
!  do
!    read(1022,*,iostat=io)
!    if (io < 0) exit
!    np = np + 1
!  enddo
!  rewind(1022)
!
!  NumForceNodes = np
!  print *,'Number of nodes with external point load BC is:',NumForceNodes
!  allocate(NodalLoads(NumNodes,6)); NodalLoads = 0.0
!  allocate(RHS(NumNodes*6,1)); RHS = 0.0
!  allocate(NodeNumb(NumForceNodes));
 
!NodalLoads(node number,1:6) is a two-dimensional array of forces and moments at all nodes in the global system
!Each row has seven columns. First col is node number, 2nd through 7th will give Fx,Fy,Fz,Mx,My,Mz in the inertial frame
! do num=1,NumForceNodes
!     read(1022,*)NodeNumb(num),NodalLoads(NodeNumb(num),1:6)
!  enddo
!  rewind(1022)
!  close(1022)
! Updating the RHS vector
!  j = 0
!  do num=1,NumNodes   !Node number loop
!        do i = 1,6    !DOF number loop
!           j = j + 1
!           RHS(j,1) = NodalLoads(num,i)
!        enddo
!  enddo
!  allocate(RHS(NumNodes*6,1)); RHS = 0.0;
!  NumForceNodes = 0
!-------------------------------------------------- 
! Read boundary condition information
!  open(1032,file='bc.dat',action='read')
!
!  np=0
!  do
!     read(1032,*,iostat=io)
!     if (io < 0) exit
!     np = np + 1
!  enddo
!
!  rewind(1032)
!  NumBC = np
!  print *,'Number of nodes with prescribed displacement/rotation B.C. is:',NumBC

!  allocate(bcdof(NumBC)); bcdof = 0
!  allocate(bcval(NumBC)); bcval = 0.0
!NodalLoads(node number,1:6) is a two-dimensional array of forces and moments at all nodes in the global system
!  do num=1,NumBC
!     read(1032,*)bcdof(num),bcval(num)
!  enddo
!  close(1032)
!  deallocate(NodalLoads)
!  deallocate(NodeNumb)
!====================================================
! For Coupled code of Loci-STREAM:

    allocate(RHS(NumNodes*6,1)); RHS = 0.0;
    NumForceNodes = 0
    NumBc = LociSTREAMNumBC;
    allocate(bcdof(NumBC));
    allocate(bcval(NumBC));
    do i=1,NumBc
      bcdof(i) = LociSTREAMBCdof(i)
      bcval(i) = LociSTREAMBCval(i)
    enddo

 end subroutine nlams_nodal_loads
end module nlams_loads 
