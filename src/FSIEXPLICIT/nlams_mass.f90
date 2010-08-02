!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_MASS. Satish Chimakurthi. 2Jan2009
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!     Compute element and global mass mass, damping, force vectors
!  
!-> Reference:
!     Chimakurthi SK, Stanford BK, Cesnik CES, Shyy W.
!     Flapping Wing CFD/CSD Aeroelastic Formulation Based on Co-rotational Shell
!     Finite Element
!     Presented in AIAA SDM 2009, AIAA 2009-2412
!
!->Subroutines:
!
!     |-nlams_massmat
!
!
!->Remarks:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nlams_mass
    
contains
  subroutine nlams_massmat
    use mesh_connect_vars,only:Xgl,Ygl,Zgl,NumElems,Connect,NumNodesPerEl,NumDofPerNode
    use nlams_undeformed
    use nlams_input,only:RHO,Thickness,alfab,TStep
    use nlams_init
    use nlams_transforms
    use lib_fem
    
    implicit none
!Local variables
    integer:: iElem,i,iGauss,j,k,p,q
    real(8):: factor,ElemArea,WeightM
    real(8),dimension(3,18) :: N    
    real(8),dimension(3,1) ::Xg,Yg,Zg,Vec1,Vec2,Vec3,Aux1,Aux2,Aux3
    real(8)::Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3       
    real(8)::xij,yij,lij,x1,x2,x3,y1,y2,y3,z1,z2,z3,y31,y12,x31,x12
    real(8)::y13,x13,y23,y21,x23,x21,y32,x32
    real(8),dimension(3) :: Gauss

!Gauss points and weights for mass matrix
    real(8) ::Coords(7,2)
    real(8) ::Weight(7)  
    real(8),save,dimension(3,1)::MidCoord,MidCoordLocal
    real(8),dimension(18,18)::Elem_LocalCSYS_Full_Init  !This matrix is simply formed by using Elem_LocalCSYS_Init along the diagonal
    real(8),save,dimension(18,18) :: M_elem,K_elem_cg 
    real(8),save,dimension(18,1)  :: Force_mbd_elem
    
!    real(8),dimension(3,3)::tempval
    real(8),dimension(3,3)::tempmat
    real(8),dimension(3,1)::tempval2
    real(8),dimension(3,9)::Nxyz
    real(8),dimension(9,1)::LocalCoordVec

    M_shell = 0.0d0
    K_cg    = 0.0d0
    Force_mbd = 0.0d0

    tempmat = 0.0d0
    tempval2= 0.0d0

!    print *, 'RHO in massmat:', RHO

!Looping through all the elements
    do iElem=1,NumElems
     
       M_elem = 0.0d0
       K_elem_cg = 0.0d0
       Force_mbd_elem  = 0.0d0

!Reading FE mesh coordinates of all the nodes in the current element in arrays Xglt, Yglt, Zglt
!Xglt,Yglt,Zglt contain x,y,z coordinates respectively, of all nodes belonging to element iElem
       
       do i=1,3
          Xg(i,1)=Xgl(Connect(iElem,i))
          Yg(i,1)=Ygl(Connect(iElem,i))
          Zg(i,1)=Zgl(Connect(iElem,i))
       enddo

! MidCoord is location of node 1 at each element at initial configuration with respect to global frame  
       MidCoord(1,1) = Xg(1,1)
       MidCoord(2,1) = Yg(1,1)
       MidCoord(3,1) = Zg(1,1)
      
       Vec1(1,1) = Xg(1,1)
       Vec1(2,1) = Yg(1,1)
       Vec1(3,1) = Zg(1,1)
       
       Vec2(1,1) = Xg(2,1)
       Vec2(2,1) = Yg(2,1)
       Vec2(3,1) = Zg(2,1)
       
       Vec3(1,1) = Xg(3,1)
       Vec3(2,1) = Yg(3,1)    
       Vec3(3,1) = Zg(3,1)
       
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

! Elem_LocalCSYS_Init will transform the component of vector in local to global frame
! Transpose of Elem_LocalCSYS_Init will transform the componet of a vector in global to local frame
       call shell_transf(Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3,Elem_LocalCSYS_Init(iElem,:,:))

! Positional vector of each node from node 1 with respect to local frame at initial time
       Vec1 = matmul(transpose(Elem_LocalCSYS_Init(iElem,:,:)),Aux1)
       Vec2 = matmul(transpose(Elem_LocalCSYS_Init(iElem,:,:)),Aux2)
       Vec3 = matmul(transpose(Elem_LocalCSYS_Init(iElem,:,:)),Aux3)
       
!x1,x2,x3 are the global x-coordinates of nodes of the current element expressed in the local element frame
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
       y23 = y2-y3
       y21 = -y12
       x23 = x2-x3
       x21 = -x12
       y32 = y3-y2
       x32 = x3-x2
       
! Area of the current element computed based on initial configuration     
       ElemArea = 0.5d0*(y12*x31 - x12*y31)
              
! Gauss points for consistent mass matrix calculation
       call fem_gauss_val(0,7,Coords,Weight)      

       do iGauss=1,7 
          Gauss(2)=          Coords(iGauss,1)
          Gauss(3)=          Coords(iGauss,2)
          Gauss(1)=1.d0 - Gauss(2) - Gauss(3)
           
          WeightM = Weight(iGauss)

          N_dkt(1,1) = ((Gauss(1)**2)*(3 - 2*Gauss(1))) + (2*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,2) = (-Gauss(1)**2)*(y12*Gauss(2) + y13*Gauss(3)) - ((0.5)*(y12 + y13)*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,3) = (Gauss(1)**2)*(x12*Gauss(2) + x13*Gauss(3)) + ((0.5)*(x12 + x13)*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,4) = ((Gauss(2)**2)*(3 - 2*Gauss(2))) + (2*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,5) = (-Gauss(2)**2)*(y23*Gauss(3) + y21*Gauss(1)) - ((0.5)*(y23 + y21)*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,6) = (Gauss(2)**2)*(x23*Gauss(3) + x21*Gauss(1)) + ((0.5)*(x23 + x21)*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,7) = ((Gauss(3)**2)*(3 - 2*Gauss(3))) + (2*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,8) = (-Gauss(3)**2)*(y31*Gauss(1) + y32*Gauss(2)) - ((0.5)*(y31 + y32)*Gauss(1)*Gauss(2)*Gauss(3))
          N_dkt(1,9) = (Gauss(3)**2)*(x31*Gauss(1) + x32*Gauss(2)) + ((0.5)*(x31 + x32)*Gauss(1)*Gauss(2)*Gauss(3))

          N_opt(1,1) = Gauss(1)
          N_opt(1,2) = 0.0d0
          N_opt(1,3) = 0.5d0*alfab*Gauss(1)*(y12*Gauss(2) - y31*Gauss(3))
          N_opt(1,4) = Gauss(2)
          N_opt(1,5) = 0.0d0
          N_opt(1,6) = 0.5d0*(alfab*Gauss(2))*(y23*Gauss(3) - y12*Gauss(1))
          N_opt(1,7) = Gauss(3)
          N_opt(1,8) = 0.0d0
          N_opt(1,9) = 0.5d0*(alfab*Gauss(3))*(y31*Gauss(1) - y23*Gauss(2))
          
          N_opt(2,1) = 0.0d0
          N_opt(2,2) = Gauss(1)
          N_opt(2,3) = 0.5d0*(alfab*Gauss(1))*(x21*Gauss(2) - x13*Gauss(3))
          N_opt(2,4) = 0.0d0
          N_opt(2,5) = Gauss(2)
          N_opt(2,6) = 0.5d0*(alfab*Gauss(2))*(x32*Gauss(3) - x21*Gauss(1))
          N_opt(2,7) = 0.0d0
          N_opt(2,8) = Gauss(3)
          N_opt(2,9) = 0.5d0*(alfab*Gauss(3))*(x13*Gauss(1) - x32*Gauss(2))

!Form the N shape function matrix
!          N = 0.0d0
! N(3,18) 
          N(1,1) = N_opt(1,1)
          N(1,2:5) = 0.0
          N(1,6) = N_opt(1,3)
          N(1,7) = N_opt(1,4)
          N(1,8) = 0.0
          N(1,9:11) = 0.0
          N(1,12) = N_opt(1,6)
          N(1,13) = N_opt(1,7)
          N(1,14) = 0.0
          N(1,15:17) = 0.0
          N(1,18) = N_opt(1,9)
 
          N(2,1) = 0.0
          N(2,2) = N_opt(2,2)
          N(2,3:5) = 0.0
          N(2,6) = N_opt(2,3)
          N(2,7) = 0.0
          N(2,8) = N_opt(2,5)
          N(2,9:11) = 0.0
          N(2,12) = N_opt(2,6)
          N(2,13) = 0.0
          N(2,14) = N_opt(2,8)
          N(2,15:17) = 0.0
          N(2,18) = N_opt(2,9)
          !
          N(3,1:2) = 0.0
          N(3,3:5) = N_dkt(1,1:3)
          N(3,6:8) = 0.0
          N(3,9:11) = N_dkt(1,4:6)
          N(3,12:14) = 0.0
          N(3,15:17) = N_dkt(1,7:9)
          N(3,18) = 0.0
            
! Eq. (1.22) in the Satish's SDM 2009 paper
!          factor = 0.0d0
! // Need to be modified for composite model: i.e. Thickness will be varied. (Thicknes(ielem))
          factor = RHO*ElemArea*Thickness

! M_elem(18,18): element local mass matrix
          M_elem = M_elem + factor*matmul(transpose(N),N)*WeightM

 
! Element gyroscopic damping matrix: Eq. (1.23) in the Satish's SDM 2009 paper
        K_elem_cg = K_elem_cg + 2.d0*factor*(WeightM*matmul(matmul(matmul(matmul(matmul(matmul(transpose(N), &
& transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),SkewOmega),T),Elem_LocalCSYS_Init(iElem,:,:)),N))
 
!        tempval =0.0
!        tempval = (SkewOmegadot + matmul(SkewOmega,SkewOmega))
        tempmat = (SkewOmegadot + matmul(SkewOmega,SkewOmega))

        Nxyz = 0.0
        Nxyz(1,1) = Gauss(1)
        Nxyz(2,2) = Gauss(1)
        Nxyz(3,3) = Gauss(1)
        Nxyz(1,4) = Gauss(2)
        Nxyz(2,5) = Gauss(2)
        Nxyz(3,6) = Gauss(2)
        Nxyz(1,7) = Gauss(3)
        Nxyz(2,8) = Gauss(3)
        Nxyz(3,9) = Gauss(3)

        LocalCoordVec(1,1) = x1
        LocalCoordVec(2,1) = y1
        LocalCoordVec(3,1) = z1
        LocalCoordVec(4,1) = x2
        LocalCoordVec(5,1) = y2
        LocalCoordVec(6,1) = z2
        LocalCoordVec(7,1) = x3
        LocalCoordVec(8,1) = y3
        LocalCoordVec(9,1) = z3

        tempval2 = matmul(matmul(Elem_LocalCSYS_Init(iElem,:,:),Nxyz),LocalCoordVec)

! Eq. (1.25) in Satish's SDM 2009 paper: The force vector is due to the prescribed rigid body motion
  !  Force_mbd_elem = Force_mbd_elem + WeightM*(factor*matmul(matmul(matmul(matmul(matmul(transpose(N),transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),tempval),T),MidCoord)) +(SkewOmegadot + matmul(SkewOmega,SkewOmega)) (WeightM*factor*matmul(matmul(matmul(transpose(N),transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),XRdot)) + (WeightM*factor*matmul(matmul(matmul(matmul(matmul(transpose(N),transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),tempval),T),tempval2))
!        Force_mbd_elem = Force_mbd_elem + term1 + term2 + term3

  Force_mbd_elem = Force_mbd_elem + WeightM*factor*matmul(matmul(matmul(matmul(matmul(transpose(N),&
&transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),tempmat),T),MidCoord) &
&                                 + WeightM*factor*matmul(matmul(matmul(transpose(N),&
&transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),XRdot)                              &
&                                 + WeightM*factor*matmul(matmul(matmul(matmul(matmul(transpose(N),&
&transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),tempmat),T),tempval2)


       enddo ! igauss

! Form Elem_LocalCSYS_Full
    k = 1
    do i=1,6
       Elem_LocalCSYS_Full_Init(k:k+2,k:k+2) = Elem_LocalCSYS_Init(iElem,1:3,1:3)
       k = k + 3
    enddo

! Transform elment matrices and vectors from local frame to the global frame
       M_elem = matmul(matmul(Elem_LocalCSYS_Full_Init,M_elem),transpose(Elem_LocalCSYS_Full_Init))
       K_elem_cg = matmul(matmul((Elem_LocalCSYS_Full_Init),K_elem_cg),transpose(Elem_LocalCSYS_Full_Init))
       Force_mbd_elem = matmul(Elem_LocalCSYS_Full_Init,-Force_mbd_elem)
 
       call nlams_index(iElem)

! Assembling 
       do i=1,NumNodesPerEl*NumDofPerNode
          p = ElDofIndex(i)
          do j=1,NumNodesPerEl*NumDofPerNode
             q = ElDofIndex(j)
             M_shell(p,q) = M_shell(p,q) + M_elem(i,j)
             K_cg(p,q) = K_cg(p,q) + K_elem_cg(i,j)
          enddo
       enddo
       
       do i=1,NumNodesPerEl*NumDofPerNode
          p = ElDofIndex(i)
          Force_mbd(p,1) = Force_mbd(p,1) + Force_mbd_elem(i,1)
       enddo
       
    enddo

!
!  open(10,file='M_shell_new.dat')
!  do i=1,AnsSize
!   write(10,*)  M_shell(i,1),  M_shell(i,2),  M_shell(i,3)
!  enddo
!  close(10)
!
!  open(10,file='K_cg_new.dat')
!  do i=1,AnsSize
!   write(10,*) K_cg(i,1),K_cg(i,2),K_cg(i,3)
!  enddo
!  close(10)
!
!  open(10,file='Force_mbd_new.dat')
!  do i=1,AnsSize
!   write(10,*) Force_mbd(i,1)
!  enddo
!  close(10)
!
!
!  if(TStep.eq.2)  stop




  end subroutine nlams_massmat
end module nlams_mass

