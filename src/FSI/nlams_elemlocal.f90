!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- ElemLocal Satish Chimakurthi. 2Jan2009
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!     Compute nodal displacment and rotation  
!
!-> Reference:
!     Chimakurthi SK, Stanford BK, Cesnik CES, Shyy W.
!     Flapping Wing CFD/CSD Aeroelastic Formulation Based on Co-rotational Shell
!     Finite Element
!     Presented in AIAA SDM 2009, AIAA 2009-2412
!
!->Subroutines:
!
!     |-Elem_LocalCoord
!
!
!->Remarks:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ElemLocal
! Global vaiables
  real(8),save::x1,x2,x3,y1,y2,y3,z1,z2,z3,y31,y12,x31,x12
  real(8),save::y13,x13,y23,y21,x23,x21,y32,x32
  real(8),save::ElemArea
  
contains
  
subroutine Elem_LocalCoord(iElem,Xg,Yg,Zg)
    use nlams_transforms
    use nlams_init, only: Elem_LocalCSYS,Nodal_CSYS,Elem_LocalCSYS_Init,ElePureD
    use mesh_connect_vars,only: NumNodesPerEl,NumElems,Xgl,Ygl,Zgl,Connect
    use interface_lapack
    use nlams_input, only: TStep, Iter
    
    implicit none
!intent(in) variables
    real(8),dimension(3,1),intent(in)::Xg,Yg,Zg
    integer,intent(in)  :: iElem
    
!Local variables
    integer:: i,j
    integer:: N, error
    real(8),dimension(3,1)::Vec1,Vec2,Vec3,Aux1,Aux2,Aux3,XgI,YgI,ZgI
    real(8)::Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3
    real(8)::xi1,xi2,xi3,yi1,yi2,yi3,zi1,zi2,zi3,u11,u12,u13,u21,u22,u23,u31,u32,u33
    real(8),dimension(3,3) :: Spin_Tens,Iden,T_Eo,T_S,T,T_e
    real(8),dimension(3,3,1) :: DefRot
    real(8),dimension(3,3) :: temp2,temp1,invtemp2

! Vec1=0.0d0
! Vec2=0.0d0
! Vec3=0.0d0
! Aux1=0.0d0
! Aux2=0.0d0
! Aux3=0.0d0

! Positional vectors at each node with respect to global frame
! Node 1
    Vec1(1,1) = Xg(1,1)
    Vec1(2,1) = Yg(1,1)
    Vec1(3,1) = Zg(1,1)
! Node 2    
    Vec2(1,1) = Xg(2,1)
    Vec2(2,1) = Yg(2,1)
    Vec2(3,1) = Zg(2,1)
! Node 3    
    Vec3(1,1) = Xg(3,1)
    Vec3(2,1) = Yg(3,1)    
    Vec3(3,1) = Zg(3,1)

! Positional vectors at each node from node 1 with respect to global frame     
    Aux1 = Vec1 - Vec1
    Aux2 = Vec2 - Vec1
    Aux3 = Vec3 - Vec1
    
    Xg1=Xg(1,1)
    Xg2=Xg(2,1)
    Xg3=Xg(3,1)
    
    Yg1=Yg(1,1)
    Yg2=Yg(2,1)
    Yg3=Yg(3,1)
    
    Zg1=Zg(1,1)
    Zg2=Zg(2,1)
    Zg3=Zg(3,1)
      
! Elem_LocalCSYS will transform the components of a vector in local frame to global frame

    call shell_transf(Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3,Elem_LocalCSYS(iElem,:,:))

! Transpose(Elem_LocalCSYS) will transoform the comonents of a vector from global frame to local frame
! Positional vectors at each nodes from node 1 with respect to local frame   
    Vec1 = matmul(transpose(Elem_LocalCSYS(iElem,:,:)),Aux1)
    Vec2 = matmul(transpose(Elem_LocalCSYS(iElem,:,:)),Aux2)
    Vec3 = matmul(transpose(Elem_LocalCSYS(iElem,:,:)),Aux3)
    
! x1,y1,z1 are coordinate system of node 1 from node 1 with respect to local frame
! x2,y2,z2 are coordinate system of node 2 from node 1 with respect to local frame
! x3,y3,z3 are coordinate system of node 3 from node 1 with respect to local frame

    x1 = Vec1(1,1)
    x2 = Vec2(1,1)
    x3 = Vec3(1,1)

    y1 = Vec1(2,1)
    y2 = Vec2(2,1)
    y3 = Vec3(2,1)

    z1 = Vec1(3,1)
    z2 = Vec2(3,1)
    z3 = Vec3(3,1)
    
    y31 = y3 - y1
    y12 = y1 - y2
    x31 = x3 - x1
    x12 = x1 - x2
    
    y13 = -y31
    x13 = -x31
    y23 = y2 - y3
    y21 = -y12
    x23 = x2 - x3
    x21 = -x12
    y32 = y3 - y2
    x32 = x3 - x2
    
!Area of the current element     
    ElemArea = 0.5d0*(y12*x31 - x12*y31)
    
!Compute pure element deformational translations and rotations reading the initial coordinates of the current element in the global frame
    do i =1,3
       XgI(i,1)=Xgl(Connect(iElem,i))
       YgI(i,1)=Ygl(Connect(iElem,i))
       ZgI(i,1)=Zgl(Connect(iElem,i))
    enddo
    
    Vec1(1,1) = XgI(1,1)
    Vec1(2,1) = YgI(1,1)
    Vec1(3,1) = ZgI(1,1)
    
    Vec2(1,1) = XgI(2,1)
    Vec2(2,1) = YgI(2,1)
    Vec2(3,1) = ZgI(2,1)
    
    Vec3(1,1) = XgI(3,1)
    Vec3(2,1) = YgI(3,1)
    Vec3(3,1) = ZgI(3,1)
    
    Aux1 = Vec1 - Vec1
    Aux2 = Vec2 - Vec1
    Aux3 = Vec3 - Vec1

    Xg1=XgI(1,1)
    Xg2=XgI(2,1)
    Xg3=XgI(3,1)

    Yg1=YgI(1,1)
    Yg2=YgI(2,1)
    Yg3=YgI(3,1)

    Zg1=ZgI(1,1)
    Zg2=ZgI(2,1)
    Zg3=ZgI(3,1)

    call shell_transf(Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3,Elem_LocalCSYS_Init(iElem,:,:))

    Vec1 = matmul(transpose(Elem_LocalCSYS_Init(iElem,:,:)),Aux1)
    Vec2 = matmul(transpose(Elem_LocalCSYS_Init(iElem,:,:)),Aux2)
    Vec3 = matmul(transpose(Elem_LocalCSYS_Init(iElem,:,:)),Aux3)

!xi1,xi2,xi3 are the global x-coordinates of nodes of the initial element expressed in the local element frame
    xi1 = Vec1(1,1)
    xi2 = Vec2(1,1)
    xi3 = Vec3(1,1)

!yi1,yi2,yi3 are the global y-coordinates of nodes of the initial element expressed in the local element frame
    yi1 = Vec1(2,1)
    yi2 = Vec2(2,1)
    yi3 = Vec3(2,1)

!zi1,zi2,zi3 are the global z-coordinates of nodes of the initial element expressed in the local element frame
    zi1 = Vec1(3,1)
    zi2 = Vec2(3,1)
    zi3 = Vec3(3,1)

!Pure element deformation translations
!Node 1
    u11 = x1 - xi1
    u12 = y1 - yi1
    u13 = z1 - zi1
    
!Node 2
    u21 = x2 - xi2
    u22 = y2 - yi2
    u23 = z2 - zi2
    
!Node 3
    u31 = x3 - xi3
    u32 = y3 - yi3
    u33 = z3 - zi3
    
!Compute: Pure element deformation rotations in local coodinate system E
    T_S  = 0.0d0
    T_Eo = 0.0d0
    T_e  = 0.0d0
    temp1= 0.0d0
    temp2= 0.0d0
    invtemp2 = 0.0d0

    do i=1,NumNodesPerEl
! Transformation matrix for trid S
       T_S  = Nodal_CSYS(i,iElem,:,:)
! Transoformation matrix for coodinate system Eo     
       T_Eo = Elem_LocalCSYS_Init(iElem,:,:)  
! Transoformation matrix for coodinate system E
       T_e  = Elem_LocalCSYS(iElem,:,:)  

!       if(TStep.gt.1.and.iElem.eq.5)then
!       do j=1,3
!        write(*,*) T_S(j,1:3)
!       enddo
!       do j=1,3
!        write(*,*) T_Eo(j,1:3)
!       enddo
!       do j=1,3
!        write(*,*) T_e(j,1:3)
!       enddo
!       read(*,*)
!       endif


!Identity matrix
       Iden = 0.d0
       Iden(1,1) = 1.0
       Iden(2,2) = 1.0
       Iden(3,3) = 1.0
       
! Eq. (46) in Khosravi's paper: An orthogonal transformation matrix which describes rotation of nodal triad So to S in local cooridnate E
! This T is local variables for computing Spin_Tens
! This is not transformation matrix T_IG which will transfrom from local to global
       T = matmul(matmul(transpose(T_e),T_S),T_Eo)
       temp1 = T - Iden
       temp2 = T + Iden
       N = 3
       call lapack_inv(N,temp2,invtemp2,error)

! Eq. (48) in Khosravi's paper
       Spin_Tens = 0.0d0
       Spin_Tens = 2.d0*matmul(temp1,invtemp2)    

! Pure nodal rotation expressed in E are equal to the components of an antisymmetric matrix Spin_Tens        
       DefRot(i,1,1) = Spin_Tens(3,2)
       DefRot(i,2,1) = Spin_Tens(1,3)
       DefRot(i,3,1) = Spin_Tens(2,1)
    enddo

!Vector of pure deformational element deformations in local coordinate system E
! The pure nodal translation in E3 direction are zero, since both corotated and 
! current configurations are coplanar.

!Node 1 of the current element
!ElePureD is d_i i.e. Equation 50 of the reference paper
    ElePureD(1,1) = u11 
    ElePureD(2,1) = u12
    ElePureD(3,1) = u13
    ElePureD(4,1) = DefRot(1,1,1) 
    ElePureD(5,1) = DefRot(1,2,1) 
    ElePureD(6,1) = DefRot(1,3,1) 
    
!Node 2 of the current element
    ElePureD(7,1)  = u21 
    ElePureD(8,1)  = u22 
    ElePureD(9,1)  = u23 
    ElePureD(10,1) = DefRot(2,1,1) 
    ElePureD(11,1) = DefRot(2,2,1) 
    ElePureD(12,1) = DefRot(2,3,1) 
    
!Node 3 of the current element
    ElePureD(13,1) = u31 
    ElePureD(14,1) = u32 
    ElePureD(15,1) = u33 
    ElePureD(16,1) = DefRot(3,1,1) 
    ElePureD(17,1) = DefRot(3,2,1) 
    ElePureD(18,1) = DefRot(3,3,1) 


!    write(23033,*)iElem,sum(ElePureD)

  end subroutine Elem_LocalCoord
  
end module ElemLocal
