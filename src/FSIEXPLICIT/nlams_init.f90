module nlams_init
 
  real(8),save,dimension(:,:,:),allocatable::Elem_LocalCSYS      !Element local coordinate systems
  real(8),save,dimension(:,:,:),allocatable::Elem_LocalCSYS_Init !Element local coordinate systems
  real(8),save,dimension(:,:,:,:),allocatable::Nodal_CSYS        !Nodal coordinate systems
  real(8),save,dimension(:,:,:,:),allocatable::Nodal_CSYS_pre    !Nodal coordinate systems

  real(8),dimension(:),allocatable::WeightGauss                  !Gauss weights for numerical integration
  real(8),dimension(:,:),allocatable::CoordGauss                 !Gauss point coordinates
  real(8),dimension(:,:),allocatable::K_shell,K_shell_n          !Global shell stiffness matrix
  real(8),dimension(:,:),allocatable::C_shell                    !Global shell damping matrix
  real(8),dimension(:,:),allocatable::M_shell                    !Global shell mass matrix
  real(8),dimension(:,:),allocatable::K_dyn,K_cg                 !Global shell dynamic damping and stiffen matrices
  real(8),dimension(:,:),allocatable::Force_mbd,Force_mbd_pre    !Global multi-body dynamic inertial force 
!  real(8),dimension(:,:),allocatable::InvM_shell                 !Global inverse shell mass matrix
  real(8),dimension(:,:),allocatable::IntForce_g,IntForce_g_n    !Global internal force vector
 
  real(8),save,allocatable:: EleStVec(:,:)                       !Element state vector
  real(8),save,allocatable:: EleStVecCurr(:,:)                   !Element state vector
  integer,save,allocatable:: ElDofIndex(:)                       !Indices
!  real(8),save,allocatable:: NodeLocalCoord(:,:,:)
  real(8),save,allocatable:: FinalAns(:)
  real(8),save,allocatable:: PreviousAns(:)
  integer  ::AnsSize

  real(8),save,dimension(:,:),allocatable :: D_pre, Ddot_pre, Dddot_pre
  real(8),save,dimension(:,:),allocatable :: D_curr,Ddot_curr,Dddot_curr,Rcorr,EqError
  real(8),save,dimension(2,9)::N_opt,dN_opt_g2,dN_opt_g3
  real(8),save,dimension(1,9)::N_dkt,dN_dkt_g2,dN_dkt_g3
  real(8),dimension(18,1) ::ElePureD

  real(8),save,allocatable :: AerForce(:,:),AerForce_pre(:,:)
  real(8),save,allocatable :: DefGrid(:,:)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !-> Subroutine NLAMS__INIT
  !
  !-> Description:
  !
  !    Allocate memory for the description of the undeformed (rigid body)
  !    motion of the system. Then, initialize the variables.
  !
  !-> Remarks.-
  !
  !  1) The undeformed location at center points of member elements is given
  !     with respect to the a-frame of the member. It has two components, 
  !     the position of the element along the X axis of the a-frame, plus
  !     the effect of the member curvature.
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine nlams_allocate_init(Num1,Num2)
    use LociSTREAMdata
    use nlams_input
    use mesh_connect_vars, only: NumElems,NumNodes,NumElDof,NumDofPerNode,NumNodesPerEl

    implicit none
    integer,intent(in)::Num1,Num2  !Num1 is NumNodesPerEl and Num2 is NumDofPerNode

!    print *,'Maximum number of inner iterations is:',MaxIter
    
! Four dimensional array containing the local element coordinate systems
    allocate(Elem_LocalCSYS(NumElems,3,3)); Elem_LocalCSYS = 0.0
    allocate(Elem_LocalCSYS_Init(NumElems,3,3)); Elem_LocalCSYS_Init = 0.0

! Nodal coordinate systems
    allocate(Nodal_CSYS(NumNodesPerEl,NumElems,3,3)); Nodal_CSYS=0.0
    allocate(Nodal_CSYS_pre(NumNodesPerEl,NumElems,3,3)); Nodal_CSYS_pre=0.0
  
! Assuming the global coordinate system is defined for only one wing with X in the spanwise diren
! Y in the chordwise direction, Z vertical, Nodal_CSYS is initialized as follows
if(TStep==1)then
     Nodal_CSYS(:,:,1,1) = 1.0
     Nodal_CSYS(:,:,1,2) = 0.0
     Nodal_CSYS(:,:,1,3) = 0.0
     Nodal_CSYS(:,:,2,1) = 0.0
     Nodal_CSYS(:,:,2,2) = 1.0
     Nodal_CSYS(:,:,2,3) = 0.0
     Nodal_CSYS(:,:,3,1) = 0.0
     Nodal_CSYS(:,:,3,2) = 0.0
     Nodal_CSYS(:,:,3,3) = 1.0
     
     Nodal_CSYS_pre(:,:,1,1) = 1.0
     Nodal_CSYS_pre(:,:,1,2) = 0.0
     Nodal_CSYS_pre(:,:,1,3) = 0.0
     Nodal_CSYS_pre(:,:,2,1) = 0.0
     Nodal_CSYS_pre(:,:,2,2) = 1.0
     Nodal_CSYS_pre(:,:,2,3) = 0.0
     Nodal_CSYS_pre(:,:,3,1) = 0.0
     Nodal_CSYS_pre(:,:,3,2) = 0.0
     Nodal_CSYS_pre(:,:,3,3) = 1.0
   else
     Nodal_CSYS = LociSTREAMNodal_CSYS
     Nodal_CSYS_pre = LociSTREAMNodal_CSYS
   endif
   
!    end if


!//     
!     Nodal_CSYS = LociSTREAMNodal_CSYS
!    if(TStep.eq.1)then
!    		Nodal_CSYS_pre = Nodal_CSYS
!    else
!        Nodal_CSYS_pre = LociSTREAMNodal_CSYS
!    endif

    AnsSize = NumNodes*NumDofPerNode
!    print *, 'AnsSize:' ,AnsSize

    ! Allocating memory of the global DKT stiffness matrix
    !Allocate memory for strain energy array
!    allocate(StrainEn(NumElems)); StrainEn = 0.0
    allocate(K_shell(AnsSize,AnsSize)); K_shell = 0.0
    allocate(K_shell_n(AnsSize,AnsSize)); K_shell_n = 0.0
!    allocate(K_shell_n1(AnsSize,AnsSize)); K_shell_n1 = 0.0

    ! Allocating memory for global mass matrix
    allocate(M_shell(AnsSize,AnsSize)); M_shell = 0.0
    allocate(C_shell(AnsSize,AnsSize)); C_shell = 0.0
    allocate(K_dyn(AnsSize,AnsSize)); K_dyn = 0.0
    allocate(K_cg(AnsSize,AnsSize)); K_cg = 0.0

    allocate(Force_mbd(AnsSize,1)); Force_mbd = 0.0
    allocate(Force_mbd_pre(AnsSize,1)); Force_mbd_pre = 0.0
!    allocate(Force_mbd_pre(AnsSize,1)); Force_mbd_pre = LociSTREAMforce_mbd

!    allocate(InvM_shell(AnsSize,AnsSize)); InvM_shell = 0.0
    allocate(ElDofIndex(AnsSize)); ElDofIndex = 0
    
    ! Allocating memory to the global internal force vector
    allocate(IntForce_g(AnsSize,1)); IntForce_g = 0.0
!    allocate(IntForce_g(AnsSize,1)); IntForce_g = LociSTREAMintforce_g
    allocate(IntForce_g_n(AnsSize,1)); IntForce_g_n = 0.0

! Allocating memory to the global internal force vector
!    allocate(P_n(AnsSize,1)); P_n = 0.0
!    allocate(P_n1(AnsSize,1)); P_n1 = 0.0
! Allocating memory to the global internal force vector
!    allocate(ExtForce_g(AnsSize,1)); ExtForce_g = 0.0
   
!Memory allocation for the Gauss points
    allocate(CoordGauss(NumGaussPoints,2));
    allocate(WeightGauss(NumGaussPoints));
    
!Allocating memory to element state vector
    allocate(EleStVec(NumElems,NumElDof)); EleStVec = 0.0

!Allocating memory to current element state vector
    allocate(EleStVecCurr(NumElems,NumElDof)); EleStVecCurr = 0.0

!Local coordinates of the all elements w.r.t. the origin of the initial local coordinate system
!    allocate(NodeLocalCoord(NumElems,NumNodesPerEl,3)); NodeLocalCoord = 0.0

    allocate(FinalAns(AnsSize)); FinalAns = 0.0
    allocate(PreviousAns(AnsSize)); PreviousAns = 0.0

!    allocate(D_ppre(AnsSize,1)); D_ppre = 0.0
!    allocate(D_pppre(AnsSize,1)); D_pppre = 0.0

    allocate(D_pre(AnsSize,1));     D_pre     = 0.0
    allocate(Ddot_pre(AnsSize,1));  Ddot_pre  = 0.0
    allocate(Dddot_pre(AnsSize,1)); Dddot_pre = 0.0

    allocate(D_curr(AnsSize,1)); D_curr     = LociSTREAMdisp
    allocate(Ddot_curr(AnsSize,1));  Ddot_curr  = LociSTREAMvelo
    allocate(Dddot_curr(AnsSize,1)); Dddot_curr = LociSTREAMacce

!    allocate(ICvel(AnsSize,1)); ICvel = 0.0
!    allocate(ICdisp(AnsSize,1)); ICdisp = 0.0
!    allocate(D_curr(AnsSize,1)); D_curr = 0.0
!    allocate(Ddot_curr(AnsSize,1)); Ddot_curr = 0.0
!    allocate(Dddot_curr(AnsSize,1)); Dddot_curr = 0.0

!   allocate(DeltaQ(AnsSize,1)); DeltaQ = 0.0

    allocate(Rcorr(AnsSize,1)); Rcorr = 0.0

!testing
    allocate(EqError(AnsSize,1)); EqError = 0.0

!Aerodynamic forces from CFD
    allocate(AerForce(AnsSize,1))    ; AerForce     = LociSTREAMaeroforce
    allocate(AerForce_pre(AnsSize,1)); AerForce_pre = LociSTREAMaeroforce_pre

    allocate(DefGrid(NumNodes,3)); DefGrid = 0.0

!    allocate(FsInterfaceXZone1Node(NumNodes)); FsInterfaceXZone1Node = 0.0
!    allocate(FsInterfaceYZone1Node(NumNodes)); FsInterfaceYZone1Node = 0.0
!    allocate(FsInterfaceZZone1Node(NumNodes)); FsInterfaceZZone1Node = 0.0
!    allocate(FsInterfaceXZone2Node(NumNodes)); FsInterfaceXZone2Node = 0.0
!    allocate(FsInterfaceYZone2Node(NumNodes)); FsInterfaceYZone2Node = 0.0
!    allocate(FsInterfaceZZone2Node(NumNodes)); FsInterfaceZZone2Node = 0.0

!memory for energy conservation
!  allocate(Force_int_mminus_g(AnsSize,1)); Force_int_mminus_g = 0.0  !2nd term in equation 8 of Relvas & Suleman
!  allocate(Force_int_mplus_g(AnsSize,1));  Force_int_mplus_g = 0.0
!  allocate(Glob2LocTransPre(NumElems,18,18)); Glob2LocTransPre = 0.0
!  allocate(Elem_LocalCSYS_Full_n(NumElems,18,18)); Elem_LocalCSYS_Full_n = 0.0
!  allocate(Force_mbd_mplus(AnsSize,1)); Force_mbd_mplus=0.0
!  allocate(KTmplus_g(AnsSize,AnsSize)); KTmplus_g = 0.0
!  allocate(KTmminus_g(AnsSize,AnsSize)); KTmminus_g = 0.0
!  allocate(StrainVec(AnsSize,1)); StrainVec = 0.0d0
! For eigenvalue analysis
! allocate(eigenvec(AnsSize)); eigenvec = 0.0d0
! allocate(eigenvalue(AnsSize,AnsSize)); eigenvalue = 0.0d0

  end subroutine nlams_allocate_init

  subroutine nlams_deallocate
  use mesh_connect_vars
  use nlams_loads
  use LociSTREAMdata
  
 !  if(LociSTREAMRank==0) print *, 'in side deallocation'
   deallocate(Xgl,Ygl,Zgl)
   deallocate(UndefGrid)
   deallocate(Connect)
   deallocate(RHS,bcdof,bcval)
   deallocate(Elem_LocalCSYS, Elem_LocalCSYS_Init)
   deallocate(Nodal_CSYS, Nodal_CSYS_pre)
   deallocate(K_shell,K_shell_n)
   deallocate(M_shell)
   deallocate(C_shell)
   deallocate(K_dyn)
   deallocate(K_cg)
   deallocate(Force_mbd)
   deallocate(Force_mbd_pre)
   deallocate(ElDofIndex)
   deallocate(IntForce_g)
   deallocate(IntForce_g_n)
   deallocate(CoordGauss)
   deallocate(WeightGauss)
   deallocate(EleStVec)
   deallocate(EleStVecCurr)
!   deallocate(NodeLocalCoord)
   deallocate(FinalAns)
   deallocate(PreviousAns)
   deallocate(D_pre)
   deallocate(Ddot_pre)
   deallocate(Dddot_pre)
   deallocate(D_curr)
   deallocate(Ddot_curr)
   deallocate(Dddot_curr)
   deallocate(Rcorr)
   deallocate(EqError)
   deallocate(AerForce)
   deallocate(AerForce_pre)
   deallocate(DefGrid)
!    if(LociSTREAMRank==0) print *, 'go out deallocation'
  end subroutine nlams_deallocate
end module nlams_init
