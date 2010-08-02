    subroutine readBC
    use maindata
    

    !Local variables
    integer:: io,np,num,i,j
!    real(8),save,allocatable::NodalLoads(:,:)
!    integer,save,allocatable::NodeNumb(:)
!    open(1022,file='loads.dat',action='read')
!    np=0
!    do
!     read(1022,*,iostat=io)
!     if (io < 0) exit
!     np = np + 1
!    enddo
!
!    rewind(1022)
!    mainNummForceNodes = np
!  print *,'Number of nodes with external point load BC is:',NumForceNodes
!
!  allocate(NodalLoads(NumNodes,6)); NodalLoads = 0.0
!  allocate(RHS(NumNodes*6,1)); RHS = 0.0
!  allocate(NodeNumb(NumForceNodes));
!  NodalLoads(node number,1:6) is a two-dimensional array of forces and moments at all nodes in the global system
!  Each row has seven columns. First col is node number, 2nd through 7th will give Fx,Fy,Fz,Mx,My,Mz in the inertial frame
!  do num=1,NumForceNodes
!     read(1022,*)NodeNumb(num),NodalLoads(NodeNumb(num),1:6)
!  enddo
!  rewind(1022)
!
!  close(1022)
!  !Updating the RHS vector
!  j = 0
!  do num=1,NumNodes   !Node number loop
!        do i = 1,6    !DOF number loop
!           j = j + 1
!           RHS(j,1) = NodalLoads(num,i)
!        enddo
!  enddo

! Read boundary condition information
       open(1032,file='bc.dat',action='read')

       np=0
       do
        read(1032,*,iostat=io)
        if (io < 0) exit
        np = np + 1
       enddo

       rewind(1032)
       mainNumBC = np

!  print *,'Number of nodes with prescribed displacement/rotation B.C. is:',NumBC

      allocate(mainBCdof(mainNumBC)); mainBCdof = 0
      allocate(mainBCval(mainNumBC)); mainBCval = 0.0

  !NodalLoads(node number,1:6) is a two-dimensional array of forces and moments at all nodes in the global system
  do num=1,mainNumBC
     read(1032,*)mainBCdof(num),mainBCval(num)
  enddo
  close(1032)

!  deallocate(NodalLoads)
!  deallocate(NodeNumb)
end subroutine readBC


