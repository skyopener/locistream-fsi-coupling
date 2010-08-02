!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_stiff
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!     The Optimal Membrane Element and Discrete Kirchoff Triangle Plate Bending 
!     Element Implementation
!  
!-> References:
!     Paper I: Corotational Nonlinear Analysis of Thin Plates and Shells Using a New Shell Element 
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 69, pp. 859-885, 2007
!
!     More detail of DKT bending element
!     Paper II: A Study of Three-Node Triangular Plate Bending Elements
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 15, pp. 1771-1812, 1980.
!
!     More detail of OPT membrane element
!     Paper III: A Study of Optimal Memberane triangles with Drilling Freedoms
!     Journal: Computer Methods in Applied Mechanics and Engineering, Vol. 192, pp. 2125-2168, 2003.
!
!->Subroutines:
!
!	|-nlams_stiffmat
!		|-nlams_shell_stiff_rearrange
!		|-nlams_intforce_vec
!		|-nlams_elem_frame_nodal_coord
!		|-nlams_corot_node_projector
!		|-nlams_stiff_modify
!
!->Remarks:
!    Variables are precisely as presented in the reference paper 1.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nlams_stiff
real(8),save,dimension(9,9)::K_elem_dkt,K_elem_opt
real(8),save,dimension(18,18)::K_elem_shell,K_elem_global,K_elem_dyn

!real(8),save,dimension(18,18):: KTmplus
!real(8),save,dimension(18,1) :: figmplus
!
real(8),save,dimension(18,1) :: IntForce_elem,IntForce_elem_global
real(8),save,dimension(18,18) :: K_elem_sig
real(8),save,dimension(3)::E1,E2,E3
real(8),save,dimension(18,18)::Proj

!real(8),save,dimension(18,1) :: StrainVec_elem_global

contains
  
  subroutine nlams_stiffmat
    use nlams_transforms
    use nlams_input
    use ElemLocal
    use mesh_connect_vars, only: Connect,NumElems,Xgl,Ygl,Zgl,NumNodesPerEl,NumDofPerNode
    use nlams_update
    use nlams_comp_abbd
    use nlams_form_shell
!    use nlams_init, only:ElePureD,dN_opt_g2,dN_opt_g3,dN_dkt_g2,dN_dkt_g3,Elem_LocalCSYS,K_shell,ElDofIndex,K_dyn,IntForce_g,AnsSize,N_dkt,N_opt,Elem_LocalCSYS_Init,Force_int_mplus_g,KTmplus_g,StrainVec
    use nlams_init, only:ElePureD,dN_opt_g2,dN_opt_g3,dN_dkt_g2,dN_dkt_g3,Elem_LocalCSYS,K_shell,ElDofIndex,K_dyn,IntForce_g,AnsSize,N_dkt,N_opt,Elem_LocalCSYS_Init
    use lib_fem
    use nlams_undeformed
    use interface_lapack 

    implicit none
!Local variables   
    integer::i,iElem,iGauss,pp,qq
    real(8),dimension(3,1)::Xg,Yg,Zg
    real(8)::xij,yij,lij,V,W,l13,l21,l32
    real(8),dimension(1,9)::temp1,temp2,temp3
    real(8),dimension(6)::P,q,t1,r
    real(8),dimension(9,1)::H1,H2,H3,H4
    real(8)::A,beta0,E12,E22,E11,E13,E23,E33
    real(8),dimension(9,3)::L
    real(8),dimension(3,3)::Te,Q1,Q2,Q3,Q4,Q5,Q6,Enat,ktheta1,ktheta2,ktheta3,ktheta,Btemp
    real(8),dimension(3,9)::Tthetau,B
    real(8),dimension(9,9)::km1,km2,km
    real(8),dimension(3,9)::Bb,Bm
    real(8),dimension(3,2)::CoordGauss
    real(8),dimension(3)  ::Gauss
    real(8):: E
    real(8),dimension(3,3)::EAe_ij

!Local variables for geometric stiffness computation
    integer :: ks,kend,k,j,Gaussp,m,n,ierr
    real(8),dimension(6,18)::G
    real(8),dimension(2,2)::Jac,InvJ,MembForTilda
    real(8),dimension(6,6)::Jt,MembForMat
    real(8),dimension(9,1)::Ldisp
    real(8),dimension(3,1)::MembForVec
    real(8)::factor,DetJ
    real(8)::Coords(7,2)
    real(8)::Weight(7)
    real(8),dimension(3,9)::Bm_mem
    real(8),dimension(18,18)       ::ShellEleMat,Ktemp,ShellEleMatForIntF,ShellStiffEleGlobal
    real(8),dimension(3,18)  :: N1
    real(8)::WeightM,factordyn

    real(8),dimension(2,2):: CheckMatrix

! Newton-Cotes quadrature: Gauss points are assumed apriori for numerical integration

!  Integration point 1
  CoordGauss(1,1) = 0.5d0; CoordGauss(1,2) = 0.0d0

!  Integration point 2
  CoordGauss(2,1) = 0.5d0; CoordGauss(2,2) = 0.5d0

!  Integration point 3
  CoordGauss(3,1) = 0.0d0; CoordGauss(3,2) = 0.5d0

  K_shell = 0.0d0
  IntForce_g = 0.0d0
  K_dyn = 0.0d0
  ShellEleMatForIntF = 0.0d0



  do iElem=1,NumElems
     K_elem_dkt = 0.0d0
     K_elem_opt = 0.0d0
     K_elem_sig = 0.0d0
     K_elem_dyn = 0.0d0

! Undefomed coordinate at each node per element.
        do i=1,3 
           Xg(i,1)=Xgl(Connect(iElem,i))
           Yg(i,1)=Ygl(Connect(iElem,i))
           Zg(i,1)=Zgl(Connect(iElem,i))
        enddo

! Updates global X,Y,Z coordinates of the current element
! Adding total displacement to undeformed mesh
        call nlams_globalcoord_update(iElem,Xg,Yg,Zg)  
!                  if (LociSTREAMRank==0) print *, 'after globalcoord_update '

! Updates local X,Y,Z coordinates of the current element and the local element coordinate system
! Element_LocalCSYS
        call Elem_LocalCoord(iElem,Xglt_new,Yglt_new,Zglt_new) !Compute element local coordinates

!         if (LociSTREAMRank==0) print *, 'after elem_localCoord'

! Gauss integration for DKT stiffness matrix
     do iGauss=1,3
        Gauss(2)=CoordGauss(iGauss,1)
        Gauss(3)=CoordGauss(iGauss,2)
        Gauss(1)=1-Gauss(2)-Gauss(3)

!Area of the current element     
        A = ElemArea
        
!DKT stiffness matrix
!Calculation of Pk, qk, tk, rk for k=4,5,6 for ij=23,31,12 (Eqn. 26 in the reference paper I)

! Intialization P, q, t1, and r
        P=0.0d0; q=0.0d0; t1=0.d0; r = 0.0d0

        do k=4,6
           if (k==4) then
              xij = x2 - x3
              yij = y2 - y3
              lij = sqrt((xij)**2 + (yij)**2)
           elseif (k==5) then
              xij = x3 - x1
              yij = y3 - y1
              lij = sqrt((xij)**2 + (yij)**2)
           elseif (k==6) then
              xij = x1 - x2
              yij = y1 - y2
              lij = sqrt((xij)**2 + (yij)**2)
           endif 
           P(k) = -6.0d0*xij     /(lij*lij)
           q(k) =  (3.d0*xij*yij)/(lij*lij)
           t1(k)= -6.d0*yij      /(lij*lij)
           r(k) =  (3.d0*yij*yij)/(lij*lij)
        enddo

! Eq (25) in Reference paper I               
        H1(1,1)=  P(6)*(1.d0 - 2.d0*Gauss(2)) + (P(5) - P(6))*Gauss(3)
        H1(2,1)=  q(6)*(1.d0 - 2.d0*Gauss(2)) - (q(5) + q(6))*Gauss(3)
        H1(3,1)= -4.d0 + (6.d0*(Gauss(2) + Gauss(3))) + r(6)*(1.d0 - 2.d0*Gauss(2)) - Gauss(3)*(r(5) + r(6))
        H1(4,1)= -P(6)*(1.d0 - 2.d0*Gauss(2)) + Gauss(3)*(P(4) + P(6))
        H1(5,1)=  q(6)*(1.d0 - 2.d0*Gauss(2)) - Gauss(3)*(q(6) - q(4))
        H1(6,1)= -2.d0 + 6.d0*Gauss(2) + r(6)*(1.d0 - 2.d0*Gauss(2)) + Gauss(3)*(r(4) - r(6))
        H1(7,1)= -Gauss(3)*(P(5) + P(4))
        H1(8,1)=  Gauss(3)*(q(4) - q(5))
        H1(9,1)= -Gauss(3)*(r(5) - r(4))
                
        H2(1,1)=  t1(6)*(1.d0 - 2.d0*Gauss(2)) + (t1(5) - t1(6))*Gauss(3)
        H2(2,1)=  1.d0 + r(6)*(1.d0 - 2.d0*Gauss(2)) - (r(5) + r(6))*Gauss(3)
        H2(3,1)= -q(6)*(1.d0 - 2.d0*Gauss(2)) + Gauss(3)*(q(5) + q(6))
        H2(4,1)= -t1(6)*(1.d0 - 2.d0*Gauss(2)) + Gauss(3)*(t1(4) + t1(6))
        H2(5,1)= -1.d0 + r(6)*(1.d0 - 2.d0*Gauss(2))+ Gauss(3)*(r(4) - r(6))
        H2(6,1)= -q(6)*(1.d0 - 2.d0*Gauss(2)) - Gauss(3)*(q(4) - q(6))
        H2(7,1)= -Gauss(3)*(t1(5) + t1(4))
        H2(8,1)=  Gauss(3)*(r(4) - r(5))
        H2(9,1)= -Gauss(3)*(q(4) - q(5))
        
        H3(1,1)= -P(5)*(1.d0 - 2.d0*Gauss(3)) - (P(6) - P(5))*Gauss(2)
        H3(2,1)=  q(5)*(1.d0 - 2.d0*Gauss(3)) - (q(5) + q(6))*Gauss(2)
        H3(3,1)= -4.d0 + 6.d0*(Gauss(2) + Gauss(3)) + r(5)*(1.d0 - 2.d0*Gauss(3)) - Gauss(2)*(r(5) + r(6))
        H3(4,1)=  Gauss(2)*(P(4) + P(6))
        H3(5,1)=  Gauss(2)*(q(4) - q(6))
        H3(6,1)= -Gauss(2)*(r(6) - r(4))
        H3(7,1)=  P(5)*(1.d0 - 2.d0*Gauss(3)) - Gauss(2)*(P(5) + P(4))
        H3(8,1)=  q(5)*(1.d0 - 2.d0*Gauss(3)) + Gauss(2)*(q(4) - q(5))
        H3(9,1)= -2.d0+(6.d0*Gauss(3)) + r(5)*(1.d0 - 2.d0*Gauss(3)) + Gauss(2)*(r(4) - r(5))
        
        H4(1,1)= -t1(5)*(1.d0 - 2.d0*Gauss(3)) - (t1(6) - t1(5))*Gauss(2)
        H4(2,1)=  1.d0 + r(5)*(1.d0 - 2.d0*Gauss(3)) - (r(5) + r(6))*Gauss(2)
        H4(3,1)= -q(5)*(1.d0 - 2.d0*Gauss(3)) + Gauss(2)*(q(5) + q(6))
        H4(4,1)=  Gauss(2)*(t1(4) + t1(6))
        H4(5,1)=  Gauss(2)*(r(4) - r(6))
        H4(6,1)= -Gauss(2)*(q(4) - q(6))
        H4(7,1)= t1(5)*(1.d0 - 2.d0*Gauss(3)) - Gauss(2)*(t1(5) + t1(4))
        H4(8,1)= -1.d0 + r(5)*(1.d0 - 2.d0*Gauss(3)) + Gauss(2)*(r(4) - r(5))
        H4(9,1)= -q(5)*(1.d0 - 2.d0*Gauss(3)) - Gauss(2)*(q(4) - q(5))
                
! Eq.(24) in referece paper I without mutiple by 0.5 element area
        temp1 =  y31*transpose(H1) + y12*transpose(H3)
        temp2 = -x31*transpose(H2) - x12*transpose(H4)
        temp3 = -x31*transpose(H1) - x12*transpose(H3) + y31*transpose(H2) + y12*transpose(H4)

        do k=1,9
           Bb(1,k) = temp1(1,k)
           Bb(2,k) = temp2(1,k)
           Bb(3,k) = temp3(1,k)
        enddo
        
! Eq.(24) in referece paper I
        Bb = (1.d0/(2.d0*A))*Bb

! Eq.(22) in referece paper I 
!        K_elem_dkt = K_elem_dkt + (2*ElemArea/6)*matmul(matmul(transpose(Bb),De_ij),Bb) ! Original
!        K_elem_dkt = K_elem_dkt + (2.d0*ElemArea/3.d0)*matmul(matmul(transpose(Bb),De_ij),Bb)       
        K_elem_dkt = K_elem_dkt + (2.d0*ElemArea/6.d0)*matmul(matmul(transpose(Bb),De_ij),Bb)

! OPT stiffness matrix
        if (FlagIso == 0) then
! Isotropic material 
           beta0 = 0.5*(1.d0-4.d0*NU12*NU12)
        else if(FlagIso == 1) then
          
!           E12 = Ae_ij(1,2)
!           E22 = Ae_ij(2,2)
!           E11 = Ae_ij(1,1)
!           E13 = Ae_ij(1,3)
!           E23 = Ae_ij(2,3)
!           E33 = Ae_ij(3,3)
!           W =    -6*E12**3 + 5*E22*E11**2 - (5*E22*E12**2) - E22*(75*E13**2 + 14*E13*E23 + 3*E23**2) + &
!                &  2*E12*(7*E13**2 + 46*E13*E23 + 7*E23**2) - E11*(5*E12**2 + 3*E13**2 - 6*E12*E22 - 5*E22**2 +&
!                &  14*E13*E23 + 75*E23**2) + (3*E11**2 + 82*E11*E22 + 3*E22**2 - 4*(6*E12**2 + 5*E13**2 - 6*E13*E23 + 5*E23**2))*E33 + &
!                &  4*(5*E11 - 6*E12 + 5*E22)*E33**2
! New:: Feb 22 2010 
!           E = E11*E22*E33 + E12*E23*E31 + E13*E21*E32 - E13*E22*E31 - E12*E21*E33 - E23*E32*E11
!           beta0 = max(256.d0*E/W-1.5d0,0.01)      
!-------------------------------------------------------     
            beta0 = 0.5*(1.d0-4.d0*NU12*NU12) ! temp
        endif

!Area
!        A = 0.0
        A = 0.5d0*(y21*x13 - x21*y13)
!        A = ElemArea      
        V = A*Thickness
        
!         if (LociSTREAMRank==0) print *, 'A and V: ',A, V

        
! Eq.(2) in referece paper I
        l13 = sqrt((x13)**2 + (y13)**2)
        l21 = sqrt((x21)**2 + (y21)**2)
        l32 = sqrt((x32)**2 + (y32)**2)

! Eq.(9) in reference paper I
! constant matrix over the elements: L (9,3)       
        L = 0.0

        L(1,1)= y23;   L(1,2)= 0.0d0; L(1,3)= x32
        L(2,1)= 0.0d0; L(2,2)= x32  ; L(2,3)= y23

        L(3,1)= (1.d0/6.d0)*alfab*y23*(y13 - y21)
        L(3,2)= (1.d0/6.d0)*alfab*x32*(x31 - x12)
        L(3,3)= (1.d0/3.d0)*alfab*(x31*y13 - x12*y21)

        L(4,1)= y31;   L(4,2)= 0.0d0;    L(4,3)= x13
        L(5,1)= 0.0d0; L(5,2)=  x13 ;    L(5,3)= y31

        L(6,1)= (1.d0/6.d0)*alfab*y31*(y21 - y32)
        L(6,2)= (1.d0/6.d0)*alfab*x13*(x12 - x23)
        L(6,3)= (1.d0/3.d0)*alfab*(x12*y21 - x23*y32)

        L(7,1)= y12;   L(7,2)= 0.0d0;    L(7,3)= x21
        L(8,1)= 0.0;   L(8,2)= x21  ;    L(8,3)= y12

        L(9,1)= (1.d0/6.d0)*alfab*y12*(y32 - y13)
        L(9,2)= (1.d0/6.d0)*alfab*x21*(x23 - x31)
        L(9,3)= (1.d0/3.d0)*alfab*(x23*y32 - x31*y13)

        L = L*(Thickness/2.d0)
     
! Eq.(10) in referece paper I
        Te = 0.0d0

        Te(1,1) = y23*y13*(l21**2)
        Te(1,2) = y31*y21*(l32**2)
        Te(1,3) = y12*y32*(l13**2)
        
        Te(2,1) = x23*x13*(l21**2)
        Te(2,2) = x31*x21*(l32**2)
        Te(2,3) = x12*x32*(l13**2)
        
        Te(3,1) = (y23*x31 + x32*y13)*(l21**2)
        Te(3,2) = (y31*x12 + x13*y21)*(l32**2)
        Te(3,3) = (y12*x23 + x21*y32)*(l13**2)
        
!        Te = (1/(4*A*A)) * Te
        Te = 1.d0/(4.d0*A*A) * Te
        
! Eq.(11) in referece paper I
        Tthetau = 0.0d0

        Tthetau(1,1) = x32
        Tthetau(1,2) = y32
        Tthetau(1,3) = 4.0d0*A
        Tthetau(1,4) = x13
        Tthetau(1,5) = y13
        Tthetau(1,6) = 0.0d0
        Tthetau(1,7) = x21
        Tthetau(1,8) = y21
        Tthetau(1,9) = 0.0d0
        
        Tthetau(2,1) = x32
        Tthetau(2,2) = y32
        Tthetau(2,3) = 0.0d0
        Tthetau(2,4) = x13
        Tthetau(2,5) = y13
        Tthetau(2,6) = 4.0d0*A
        Tthetau(2,7) = x21
        Tthetau(2,8) = y21
        Tthetau(2,9) = 0.0d0
        
        Tthetau(3,1) = x32
        Tthetau(3,2) = y32
        Tthetau(3,3) = 0.0d0
        Tthetau(3,4) = x13
        Tthetau(3,5) = y13
        Tthetau(3,6) = 0.0d0
        Tthetau(3,7) = x21
        Tthetau(3,8) = y21
        Tthetau(3,9) = 4.0d0*A
        
        Tthetau = (1.d0/(4.d0*A))*Tthetau

! Eq. (12) in referece paper I
        Q1 = 0.d0

        Q1(1,1) = (beta1)/(l21**2)
        Q1(1,2) = (beta2)/(l21**2)
        Q1(1,3) = (beta3)/(l21**2)
        
        Q1(2,1) = (beta4)/(l32**2)
        Q1(2,2) = (beta5)/(l32**2)
        Q1(2,3) = (beta6)/(l32**2)
        
        Q1(3,1) = (beta7)/(l13**2)
        Q1(3,2) = (beta8)/(l13**2)
        Q1(3,3) = (beta9)/(l13**2)
       
        Q1 = (2.d0*A/3.d0)*Q1
         
! Eq. (12) in referece paper I (Q2 matrix)
        Q2 = 0.d0

        Q2(1,1) = (beta9)/(l21**2)
        Q2(1,2) = (beta7)/(l21**2)
        Q2(1,3) = (beta8)/(l21**2)
        
        Q2(2,1) = (beta3)/(l32**2)
        Q2(2,2) = (beta1)/(l32**2)
        Q2(2,3) = (beta2)/(l32**2)
        
        Q2(3,1) = (beta6)/(l13**2)
        Q2(3,2) = (beta4)/(l13**2)
        Q2(3,3) = (beta5)/(l13**2)
        
        Q2 = (2.d0*A/3.d0)*Q2
      
! Eq. (12) in referece paper I (Q3 matrix)
        Q3 = 0.d0

        Q3(1,1) = (beta5)/(l21**2)
        Q3(1,2) = (beta6)/(l21**2)
        Q3(1,3) = (beta4)/(l21**2)
        
        Q3(2,1) = (beta8)/(l32**2)
        Q3(2,2) = (beta9)/(l32**2)
        Q3(2,3) = (beta7)/(l32**2)
        
        Q3(3,1) = (beta2)/(l13**2)
        Q3(3,2) = (beta3)/(l13**2)
        Q3(3,3) = (beta1)/(l13**2)
  
        Q3 = (2.d0*A/3.d0)*Q3
      
! Eq. (15) in referece paper I (Q4,Q5,Q6 matrices)
        Q4 = 0.d0
        Q5 = 0.d0
        Q6 = 0.d0

        Q4 = 0.5d0*(Q1 + Q2)
        Q5 = 0.5d0*(Q2 + Q3)
        Q6 = 0.5d0*(Q3 + Q1)
        
! Eq. (15) in referece paper I (Enat matrix)
! Ae_ij is elastic matrix not membrane stiffness matrix
        Enat = 0.0d0
        Enat = matmul(matmul(transpose(Te),Ae_ij),Te)
        
! Eq. (14) in referce paper I (ktheta matrix)
        ktheta1 = 0.d0
        ktheta2 = 0.d0
        ktheta3 = 0.d0
        ktheta  = 0.d0

        ktheta1 = matmul(matmul(transpose(Q4),Enat),Q4)
        ktheta2 = matmul(matmul(transpose(Q5),Enat),Q5)
        ktheta3 = matmul(matmul(transpose(Q6),Enat),Q6)
        ktheta  = V*(ktheta1 + ktheta2 + ktheta3)
        
! Eq. (8) in referece paper I (B matrix)
        Btemp = 0.0d0
        B     = 0.0d0
        Btemp  =  (Q1*Gauss(1) + Q2*Gauss(2) + Q3*Gauss(3))
        B  =  matmul(matmul(Te,Btemp),Tthetau)
        
! Eq. (13) in referece paper I (Km final form)
        km1 = 0.d0 ! first term
        km2 = 0.d0 ! second term
        km  = 0.d0 ! OPT stiffness matrix

! Ae_ij is elastic matrix not membrane stiffness matrix
        km1 = (matmul(matmul(L,Ae_ij),transpose(L)))/V
        km2 = (3.d0/4.d0)*beta0*matmul(matmul(transpose(Tthetau),ktheta),Tthetau)
        km = km1 + km2
        
        K_elem_opt = km
        
!Integral form of stiffness matrix as shown in Page 867 of the reference paper I (int (Bm'*E*Bm dv))
        Bm = 0.d0
        Bm =((transpose(L))/V) + (1.5d0)*sqrt(beta0)*matmul(matmul(Te,Btemp),Tthetau)

     enddo ! igauss

! This subroutine in nlams_shell.f90: K_elem_shell is 
     call nlams_shell(Bm,Bb,K_elem_opt,K_elem_dkt,K_elem_shell)

!     print *, 'After nlams_shell'
 
! In order to have the tangent stiffness matrix, it is necessary to add the stress stiffness
! (also called as geometric stiffness) to the linear conventional stiffness matrix. An expression
! for stress stiffness K_shell_sig is obtained by considering the work done by the membrane force
! Nx, Ny, Nxy as they act through displacements associated with small lateral and in-plane deflections

! Stress stiffening matrix - using 7 point Gauss integration
! StressStiff is SSTIFF of input file
     if (StressStiff == 1) then

! This subroutine in lib_fem.f90: 
! Output: Coords and Weight function for gauss integration
        call fem_gauss_val (0,7,Coords,Weight)
        
        do iGauss = 1,7 
           Gauss(2)=      Coords(iGauss,1)
           Gauss(3)=      Coords(iGauss,2)
           Gauss(1)=1.0d0-Gauss(2)-Gauss(3)
           WeightM = Weight(iGauss)

           call dkt_shape (Gauss)
           call opt_shape (Gauss)  
           
! Eq.(38) in referece paper: G matrix - 6x18 dimension
! Components of G matrix are derivates of shape function 
! with respect to g2 and g3:
           do i=1,9
              G(1,i) = dN_opt_g2(1,i) !Shape function of membrane 
              G(2,i) = dN_opt_g3(1,i) !Shape function of membrane
              G(3,i) = dN_opt_g2(2,i)
              G(4,i) = dN_opt_g3(2,i)
              G(5,i) = 0.0d0
              G(6,i) = 0.0d0
           enddo
           do i=10,18
              G(1,i) = 0.0d0
              G(2,i) = 0.0d0
              G(3,i) = 0.0d0
              G(4,i) = 0.0d0
           enddo
 
           j = 1
           do i=10,18
              G(5,i) = dN_dkt_g2(1,j)
              G(6,i) = dN_dkt_g3(1,j)
              j = j + 1
           enddo
           
! Eq. (37) in referece paper
! J matrix:
           Jac(1,1) = dN_opt_g2(1,1)*x1 + dN_opt_g2(1,4)*x2 + dN_opt_g2(1,7)*x3
           Jac(1,2) = dN_opt_g2(1,1)*y1 + dN_opt_g2(1,4)*y2 + dN_opt_g2(1,7)*y3
           Jac(2,1) = dN_opt_g3(1,1)*x1 + dN_opt_g3(1,4)*x2 + dN_opt_g3(1,7)*x3 
           Jac(2,2) = dN_opt_g3(1,1)*y1 + dN_opt_g3(1,4)*y2 + dN_opt_g3(1,7)*y3

!           DetJ = 0.d0
!           DetJ = abs(Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1))
! Inversed J matrix           
!           InvJ(1,1) =  (1.d0/DetJ)*Jac(2,2)
!           InvJ(1,2) = -(1.d0/DetJ)*Jac(1,2)
!           InvJ(2,1) = -(1.d0/DetJ)*Jac(2,1)
!           InvJ(2,2) =  (1.d0/DetJ)*Jac(1,1)
           InvJ = 0.0d0
           call lapack_inv(2,Jac,InvJ,ierr)
           
           Jt = 0.0d0
           ks = 1; kend = 2; m = 1; n = 1
           do k=1,3
              do i=ks,kend
                 do j=ks,kend
                    Jt(i,j) = InvJ(m,n)
                    n = n + 1
                 enddo
                 m = m + 1
                 n = 1
              enddo
              m = 1
              ks = ks + 2
              kend = kend + 2
           enddo
           
!           CheckMatrix = 0.0d0; 
!           CheckMatrix = matmul(InvJ,Jac)
!           write(*,*) "J and J inv"
!           do i=1,2
!           write(*,'(6e14.7)') (Jac(i,j),j=1,2),(CheckMatrix(i,j),j=1,2)
!           enddo
!           do i=1,6
!           write(*,'(6e14.7)')(Jt(i,j),j=1,6)
!           enddo
!           read(*,*)
 
! Eq.(3) in the referece paper I           
! Ldisp of size 9x1  gives only membrane pure deformation (u,v,thetaz at each of the nodes in an element)
           Ldisp = 0.0d0
           Ldisp(1:2,1) = ElePureD(1:2,1)
           Ldisp(3,1)   = ElePureD(6,1)
           Ldisp(4:5,1) = ElePureD(7:8,1)
           Ldisp(6,1)   = ElePureD(12,1)
           Ldisp(7:8,1) = ElePureD(13:14,1)
           Ldisp(9,1)   = ElePureD(18,1)
           
!Membrane strain-displacement matrix (Bm_mem) - the calculation of this was done in nlams_opt but the
!integration points used in that case were only 3. Since a 7 point integration scheme is being used
!here to evaluate the stress-stiffening matrix, Bm is being calculated again here.
!Difference between Btemp in OPT
           Btemp  = 0.d0
           Bm_mem = 0.d0

           Btemp  =  (Q1*Gauss(1) + Q2*Gauss(2) + Q3*Gauss(3))    
           Bm_mem =((transpose(L))/V) + (1.5d0)*sqrt(beta0)*matmul(matmul(Te,Btemp),Tthetau)
           
! Eq. (31) in the referece paper I
! Force resultant for membrane component   
           MembForVec = Thickness*matmul(matmul(Ae_ij,Bm_mem),Ldisp)
           
           MembForTilda(1,1) = MembForVec(1,1)
           MembForTilda(1,2) = MembForVec(3,1)
           MembForTilda(2,1) = MembForVec(3,1)
           MembForTilda(2,2) = MembForVec(2,1)
                     
! Eq. (36) in referece paper I Membrane forces matrix, Ntilda
           MembForMat = 0.0d0
           ks = 1; kend = 2; m = 1; n = 1
           do k=1,3
              do i=ks,kend
                 do j=ks,kend
                    MembForMat(i,j) = MembForTilda(m,n)
                    n = n + 1
                 enddo
                 m = m + 1
                 n = 1
              enddo
              m = 1
              ks = ks + 2
              kend = kend + 2
           enddo           
           
!           write(*,*) "MembForMat"
!           do i=1,6
!            write(*,'(6e14.7)') (MembForMat(i,j),j=1,6)
!           enddo
!           read(*,*)

! Eq.(41) in referece paper I
         factor = 2.0d0*ElemArea
         K_elem_sig = K_elem_sig +  factor*0.5*Weight(iGauss)*matmul(matmul(matmul(matmul(transpose(G),transpose(Jt)),MembForMat),Jt),G)  !think abt why 0.5 factor should be here.

!         K_elem_sig = K_elem_sig +  2.d0*A*Weight(iGauss)*matmul(matmul(matmul(matmul(transpose(G),transpose(Jt)),MembForMat),Jt),G) 

!Element dynamic stiffness matrix form the N shape function matrix
          N1 = 0.d0

          N1(1,1)    = N_opt(1,1)
          N1(1,2:5)  = 0.0d0
          N1(1,6)    = N_opt(1,3)
          N1(1,7)    = N_opt(1,4)
          N1(1,8)    = 0.0d0
          N1(1,9:11) = 0.0d0
          N1(1,12)   = N_opt(1,6)
          N1(1,13)   = N_opt(1,7)
          N1(1,14)   = 0.0d0
          N1(1,15:17)= 0.0d0
          N1(1,18)   = N_opt(1,9)

          N1(2,1)    = 0.0d0
          N1(2,2)    = N_opt(2,2)
          N1(2,3:5)  = 0.0d0
          N1(2,6)    = N_opt(2,3)
          N1(2,7)    = 0.0d0
          N1(2,8)    = N_opt(2,5)
          N1(2,9:11) = 0.0d0
          N1(2,12)   = N_opt(2,6)
          N1(2,13)   = 0.0d0
          N1(2,14)   = N_opt(2,8)
          N1(2,15:17)= 0.0d0
          N1(2,18)   = N_opt(2,9)

          N1(3,1:2)  = 0.0d0
          N1(3,3 )  = N_dkt(1,1)
          N1(3,4)   = N_dkt(1,2)
          N1(3,5)   = N_dkt(1,3)
          N1(3,6:8)  = 0.0d0
          N1(3,9:)  = N_dkt(1,4)
          N1(3,10)  = N_dkt(1,5)
          N1(3,11)  = N_dkt(1,6)
          N1(3,12:14)= 0.0d0
          N1(3,15)= N_dkt(1,7)
          N1(3,16)= N_dkt(1,8)
          N1(3,17)= N_dkt(1,9)
          N1(3,18)   = 0.0d0

          K_elem_dyn = 0.d0

!// Need to be modified: composite:
          factordyn  = RHO*ElemArea*Thickness

! Eq. (1.24) in Satish's SDM 2009 paper
          K_elem_dyn = K_elem_dyn + WeightM*(factordyn*matmul(matmul(matmul(matmul(matmul(matmul(transpose(N1),&
&transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),(SkewOmegadot+matmul(SkewOmega,SkewOmega))),T),&
&Elem_LocalCSYS_Init(iElem,:,:)),N1))
        enddo ! igauss
      endif !StressStiff == 1

! Move to outside from if(stressstiff==1)
!      do iGauss = 1,7
!           Gauss(2)=      Coords(iGauss,1)
!           Gauss(3)=      Coords(iGauss,2)
!           Gauss(1)=1.0d0-Gauss(2)-Gauss(3)
!           WeightM = Weight(iGauss)
!           call dkt_shape (Gauss)
!           call opt_shape (Gauss)
!Element dynamic stiffness matrix form the N shape function matrix
!          N1 = 0.d0
!
!          N1(1,1)    = N_opt(1,1)
!          N1(1,2:5)  = 0.0d0
!          N1(1,6)    = N_opt(1,3)
!          N1(1,7)    = N_opt(1,4)
!          N1(1,8)    = 0.0d0
!          N1(1,9:11) = 0.0d0
!          N1(1,12)   = N_opt(1,6)
!          N1(1,13)   = N_opt(1,7)
!          N1(1,14)   = 0.0d0
!          N1(1,15:17)= 0.0d0
!          N1(1,18)   = N_opt(1,9)
!
!          N1(2,1)    = 0.0d0
!          N1(2,2)    = N_opt(2,2)
!          N1(2,3:5)  = 0.0d0
!          N1(2,6)    = N_opt(2,3)
!          N1(2,7)    = 0.0d0
!          N1(2,8)    = N_opt(2,5)
!          N1(2,9:11) = 0.0d0
!          N1(2,12)   = N_opt(2,6)
!          N1(2,13)   = 0.0d0
!          N1(2,14)   = N_opt(2,8)
!          N1(2,15:17)= 0.0d0
!          N1(2,18)   = N_opt(2,9)
!
!          N1(3,1:2)  = 0.0d0
!          N1(3,3 )  = N_dkt(1,1)
!          N1(3,4)   = N_dkt(1,2)
!          N1(3,5)   = N_dkt(1,3)
!          N1(3,6:8)  = 0.0d0
!          N1(3,9:)  = N_dkt(1,4)
!          N1(3,10)  = N_dkt(1,5)
!          N1(3,11)  = N_dkt(1,6)
!          N1(3,12:14)= 0.0d0
!i          N1(3,15)= N_dkt(1,7)
!          N1(3,16)= N_dkt(1,8)
!          N1(3,17)= N_dkt(1,9)
!          N1(3,18)   = 0.0d0
!
!          K_elem_dyn = 0.d0
!
!// Need to be modified: composite:
!          factordyn  = RHO*ElemArea*Thickness
!
! Eq. (1.24) in Satish's SDM 2009 paper
!          K_elem_dyn = K_elem_dyn + WeightM*(factordyn*matmul(matmul(matmul(matmul(matmul(matmul(transpose(N1),&
!&transpose(Elem_LocalCSYS_Init(iElem,:,:))),transpose(T)),(SkewOmegadot+matmul(SkewOmega,SkewOmega))),T),&
!&Elem_LocalCSYS_Init(iElem,:,:)),N1))
!        enddo ! igauss
!     endif
!-----------------------------------------------------

! K_elem_shell= K_elem_opt + K_elem_dkt:
      
     Ktemp = K_elem_shell + K_elem_sig

! Matrix is rearranged based on the global DOF vector
! Follows Eq. (27) in referece paper I, Eq. (1.40) in Satish's SDM 2009 paper
     call nlams_shell_stiff_rearrange(Ktemp,ShellEleMat)

! Local shell stiffness matrix is rearranged based on the global DOF vector
     call nlams_shell_stiff_rearrange(K_elem_shell,ShellEleMatForIntF)

! Why ShellEleMat because this is 
! Compute the internal force vector from the rearranged local element stiffness matrix
     call nlams_intforce_vec(ShellEleMatForIntF)

! Calculate the coordinates of each node in element EleNum w.r.t. the local element frame
     call nlams_elem_frame_nodal_coord()

! Calculate the element projector matrix 
     call nlams_corot_node_projector

! Modify both stiffness and internal force vectors with the projector matrix and transform them into global frame
     call nlams_stiff_force_modify(ShellEleMat,K_elem_dyn,iElem)

! Eq. (1.43) in Satish's SDM 2009 paper
     K_elem_global = ShellEleMat
     
     call nlams_index(iElem)
!Shell stiffness matrix assembly
!ElDofIndex(:) has indices for the degree of freedom of the current element
     do i=1 ,NumNodesPerEl*NumDofPerNode
        pp = ElDofIndex(i)
        do j=1,NumNodesPerEl*NumDofPerNode
           qq = ElDofIndex(j)
           K_shell(pp,qq) = K_shell(pp,qq) + K_elem_global(i,j)
        enddo
     enddo
 
!Internal force vector assembly
     do i=1,NumNodesPerEl*NumDofPerNode
        pp = ElDofIndex(i)
        IntForce_g(pp,1) = IntForce_g (pp,1) + IntForce_elem_global(i,1)
!        StrainVec(pp,1)  = StrainVec(pp,1)    + StrainVec_elem_global(i,1)
     enddo

!For energy conservation algorithm
     !KTmplus assembly
!      do i=1,NumNodesPerEl*NumDofPerNode
!        pp = ElDofIndex(i)
!        do j=1,NumNodesPerEl*NumDofPerNode
!           qq = ElDofIndex(j)
!           KTmplus_g(pp,qq) = KTmplus_g(pp,qq) + KTmplus(i,j)
!        enddo
!     enddo
!Internal force vector assembly - for energy conservation
!     do i=1,NumNodesPerEl*NumDofPerNode
!        pp = ElDofIndex(i)
!        Force_int_mplus_g (pp,1) = Force_int_mplus_g (pp,1) + figmplus (i,1)
!     enddo

  enddo  !Elements loop
 
end subroutine nlams_stiffmat

subroutine nlams_shell_stiff_rearrange(tempin,tempout)
  !Assuming the FE is a triangle
  implicit none
  !Local variables
  real(8),dimension(18,18),intent(in) ::tempin
  real(8),dimension(18,18),intent(out)::tempout
  real(8),dimension(18,18)            ::tempo

  tempout = 0.0
  tempo   = 0.0

  tempout(1,:) = tempin(1,:)
  tempout(2,:) = tempin(2,:)
  tempout(3,:) = tempin(10,:)
  tempout(4,:) = tempin(11,:)
  tempout(5,:) = tempin(12,:)
  tempout(6,:) = tempin(3,:)
  
  tempout(7,:) = tempin(4,:)
  tempout(8,:) = tempin(5,:)
  tempout(9,:) = tempin(13,:)
  tempout(10,:) = tempin(14,:)
  tempout(11,:) = tempin(15,:)
  tempout(12,:) = tempin(6,:)
  
  tempout(13,:) = tempin(7,:)
  tempout(14,:) = tempin(8,:)
  tempout(15,:) = tempin(16,:)
  tempout(16,:) = tempin(17,:)
  tempout(17,:) = tempin(18,:)
  tempout(18,:) = tempin(9,:)
  
! Rearrange of matrix
  tempo = tempout
  tempout(:,3) = tempo(:,10)
  tempout(:,4) = tempo(:,11)
  tempout(:,5) = tempo(:,12)
  tempout(:,6) = tempo(:,3)
  tempout(:,7) = tempo(:,4)
  tempout(:,8) = tempo(:,5)
  tempout(:,9) = tempo(:,13)
  tempout(:,10) = tempo(:,14)
  tempout(:,11) = tempo(:,15)
  tempout(:,12) = tempo(:,6)
  tempout(:,13) = tempo(:,7)
  tempout(:,14) = tempo(:,8)
  tempout(:,15) = tempo(:,16)
  tempout(:,16) = tempo(:,17)
  tempout(:,17) = tempo(:,18)
  tempout(:,18) = tempo(:,9)
  
end subroutine nlams_shell_stiff_rearrange

subroutine nlams_intforce_vec(temp)
  use nlams_init, only: ElePureD
  implicit none  
  real(8),dimension(18,18),intent(in)::temp

! Eq.(1.40) in Satish's SDM 2009 paper 
 IntForce_elem = matmul(temp,ElePureD)
  
end subroutine nlams_intforce_vec

!Evaluate element nodal coordinates w.r.t. the origin of the local frame
subroutine nlams_elem_frame_nodal_coord
  use ElemLocal,only: x1,x2,x3,y1,y2,y3,z1,z2,z3
  implicit none
 
  E1 = 0.0d0; E2 = 0.0d0;  E3 = 0.0d0

  E1(1) = x1
  E1(2) = x2
  E1(3) = x3
  
  E2(1) = y1
  E2(2) = y2
  E2(3) = y3
 
  E3(1) = z1
  E3(2) = z2
  E3(3) = z3
  
end subroutine nlams_elem_frame_nodal_coord

!Calculate the projector matrix
subroutine nlams_corot_node_projector

  implicit none  
!Local variables
  integer :: i
  real(8),dimension(3,18)::Psit
  real(8),dimension(18,3)::Psi
  real(8),dimension(3,18)::Gammat
  real(8),dimension(18,3)::Gamma
  real(8),dimension(18,18)::Ide      !Identity matrix (size 18 x 18)

  Psit   = 0.0d0
  Psi    = 0.0d0
  Gammat = 0.0d0
  Gamma  = 0.0d0
  Ide    = 0.0d0

! Eq. (55) in referece paper I 
  Psit(1,1:3)   = 0.0d0
  Psit(1,4)     = 1.0d0
  Psit(1,5:9)   = 0.0d0
  Psit(1,10)    = 1.0d0
  Psit(1,11:14) = 0.0d0
  Psit(1,15)    = E2(3)
  Psit(1,16)    = 1.0d0
  Psit(1,17:18) = 0.0d0
  
  Psit(2,1:4)   = 0.0d0
  Psit(2,5)     = 1.0d0
  Psit(2,6:8)   = 0.0d0
  Psit(2,9)     = -E1(2)
  Psit(2,10)    = 0.0d0
  Psit(2,11)    = 1.0d0
  Psit(2,12:14) = 0.0d0
  Psit(2,15)    = -E1(3)
  Psit(2,16)    = 0.0d0
  Psit(2,17)    = 1.0d0
  Psit(2,18)    = 0.0d0
  
  Psit(3,1:5)   = 0.0d0
  Psit(3,6)     = 1.0d0
  Psit(3,7)     = 0.0d0
  Psit(3,8)     =  E1(2)
  Psit(3,9:11)  = 0.0d0
  Psit(3,12)    = 1.0d0
  Psit(3,13)    = -E2(3)
  Psit(3,14)    =  E1(3)
  Psit(3,15:17) = 0.0d0
  Psit(3,18)    = 1.0d0
  
  Psi = transpose(Psit)
  
  !Gamma matrix - equation 56 of Khosravi's paper
  Gammat(1,1:2)   =  0.0d0
  Gammat(1,3)     = (E1(3) - E1(2))/(E2(3)*E1(2))
  Gammat(1,4:8)   =  0.0d0
  Gammat(1,9)     = -E1(3)/(E2(3)*E1(2))
  Gammat(1,10:14) =  0.0d0
  Gammat(1,15)    =  1.0d0 /(E2(3))
  Gammat(1,16:18) =  0.0d0
  
  Gammat(2,1:2)   =  0.0d0
  Gammat(2,3)     =  1.0d0 / E1(2)
  Gammat(2,4:8)   =  0.0d0
  Gammat(2,9)     = -1.0d0 / E1(2) 
  Gammat(2,10:18) =  0.0d0
  
  Gammat(3,1)     =  0.0d0
  Gammat(3,2)     = -1.0d0 / E1(2)
  Gammat(3,4:7)   =  0.0d0
  Gammat(3,8)     =  1.0d0 / E1(2)
  Gammat(3,9:18)  =  0.0d0
  
!Ide is an 18x18 identity matrix
  do i=1,18
     Ide(i,i) = 1.0d0
  enddo
  
  Proj = 0.0d0  
! Eq.(55) in referece paper I
  Proj = Ide - matmul(Psi,Gammat)
    
end subroutine nlams_corot_node_projector

! Subroutine to project stiffness matrix and internal force vector and their transformation into global frame
! temp is ShellEleMat, tempdyn is K_elem_dyn
subroutine nlams_stiff_force_modify(temp,tempdyn,iElem)
!    use nlams_init, only: Elem_LocalCSYS,Elem_LocalCSYS_Init,ElePureD,Glob2LocTransPre,AnsSize,Elem_LocalCSYS_Full_n
    
   use nlams_init, only: Elem_LocalCSYS,Elem_LocalCSYS_Init,ElePureD,AnsSize
   implicit none
!Local variables
    integer:: i,j,k
    real(8),dimension(18,18)::K_shell_elem_project
    real(8),dimension(18,18)::Elem_LocalCSYS_Full  !This matrix is simply formed by using Elem_LocalCSYS along the diagonal
    real(8),dimension(18,1)::IntForce_elem_project,IntForce_dyn,IntForce_global_dyn
!    real(8),dimension(18,18)::Elem_LocalCSYS_Full_Init1
!    real(8),dimension(18,1)::IntForce_elem_project,IntForce_dyn,IntForce_global_dyn,IntForce_dyn_elem_project 
!    real(8),dimension(18,1)::IntForce_elem_project,IntForce_dyn,IntForce_global_dyn,IntForce_dyn_elem_project 

!intent(in) variables
    integer, intent(in):: iElem
    real(8),dimension(18,18),intent(inout)::temp
    real(8),dimension(18,18),intent(in)::tempdyn
    real(8),dimension(18,18)::tempdyn1
    
!For energy
!    real(8)::StrainEne(1,1)
!    real(8)::StrainVec_elem_global(18,1)
!    real(8),dimension(AnsSize,1)::fiLmplus
!    real(8),dimension(18,18)::Tcurr,Tpre

    tempdyn1 = 0.d0
    K_shell_elem_project      = 0.0d0
    Elem_LocalCSYS_Full       = 0.0d0
!    Elem_LocalCSYS_Full_Init1 = 0.0d0
    IntForce_elem_project     = 0.0d0
    IntForce_dyn              = 0.0d0
    IntForce_global_dyn       = 0.0d0
!    IntForce_dyn_elem_project = 0.0d0
!    StrainVec_elem_global     = 0.0d0

! Eq. (53) in referece paper I
    K_shell_elem_project   = matmul(matmul(transpose(Proj),temp),Proj)

! Eq. (52) in referece paper I    
    IntForce_elem_project  = matmul(transpose(Proj),IntForce_elem)
    
! Overwrite temp: K_shell_elem is now the projected stiffness matrix
    temp = K_shell_elem_project   
 
!Strain energy for the current element
!    StrainEne = 0.5*matmul(transpose(ElePureD),matmul(temp,ElePureD))
!    StrainEn(iElem) = StrainEne(1,1)

!Form Elem_LocalCSYS_Full
    k = 1
    do i=1,6
       Elem_LocalCSYS_Full(k:k+2,k:k+2) = Elem_LocalCSYS(iElem,1:3,1:3)
       k = k + 3
    enddo

!Compute strain distrubution on the wing
!     StrainVec_elem_global= matmul(Elem_LocalCSYS_Full,matmul(temp,ElePureD))
    
!Form Elem_LocalCSYS_Full_Init1
!    k = 1
!    do i=1,6
!       Elem_LocalCSYS_Full_Init1(k:k+2,k:k+2) = Elem_LocalCSYS_Init(iElem,1:3,1:3)
!       k = k + 3
!    enddo

! Global element stiffness matrix computation: Transform the shell element stiffness matrix into global frame
! Eq. (1.41) in Satish's SDM 2009 paper
    temp = matmul(matmul(Elem_LocalCSYS_Full,temp),transpose(Elem_LocalCSYS_Full))

! Element dynamic stiffness matrix computation in global frame
    tempdyn1 = matmul(matmul(Elem_LocalCSYS_Full,tempdyn),transpose(Elem_LocalCSYS_Full))

! Eq. (1.43) in Satish's SDM 2009 paper 
    temp = temp + tempdyn1  !Includes both shell element stiffness matrix and the dynamic stiffness term

! Transform the internal force vector into global frame
! Eq.(1.46) in Satish's SDM 2009 paper 
    IntForce_elem_global = matmul(Elem_LocalCSYS_Full,IntForce_elem_project)
 
! tempdyn is the local dynamic stiffness matrix - local dynamic internal force vector
    IntForce_dyn = matmul(tempdyn,ElePureD)

! Transform local dynamic internal force vector into global frame
    IntForce_global_dyn = matmul(Elem_LocalCSYS_Full,IntForce_dyn)

!Adding contribution of dynamic stiffness matrix
    IntForce_elem_global = IntForce_elem_global + IntForce_global_dyn
    
!-----For energy conservation algorithm based on Relvas and Suleman's 2006 reference paper--------------
!   Tcurr= transpose(Elem_LocalCSYS_Full)  !Global to local transformation matrix for the current element
!   Tpre = transpose(Elem_LocalCSYS_Full_n(iElem,:,:))
!   KTmplus = (0.4*matmul(matmul(transpose((Tpre+Tcurr)),temp),Tcurr))   !Equation 13 of Relvas and Suleman
!   fiLmplus = 0.5*(IntForce_elem_project_n + IntForce_elem_project)
!   figmplus = 0.5*matmul(transpose(Tpre+Tcurr),fiLmplus)
!----End of energy conservation changes-----------------------------------------------------------------
!  Glob2LocTransPre(iElem,:,:) = Elem_LocalCSYS_Full

 end subroutine nlams_stiff_force_modify

end module nlams_stiff


