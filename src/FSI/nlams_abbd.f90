!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_COMP_ABBD. Satish Chimakurthi. 22Aug2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   Calculation of A,B, and D matrices in [A,B;B,D]
!
!-> Reference:
!     Paper: Corotational Nonlinear Analysis of Thin Plates and Shells Using a New Shell Element
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 69, pp. 859-885, 2007
!     The above reference refers to another journal paper whose details are below:
!     Paper: A Study of Three-Node Triangular Plate Bending Elements
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 15, pp. 1771-1812, 1980.
!
!->Subroutines:
!
!     |-nlams_abbd
!
!
!->Remarks:
!    Variables are precisely as presented in the reference paper.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_comp_abbd
real(8),save,dimension(3,3)::Ae_ij,Be_ij,De_ij

contains
 subroutine nlams_abbd
!  use nlams_transforms ! ply_transf
  use nlams_input      ! E1, E2, NU12, NU21, Thickness, ply_angle, FlagIso, NumLayers

 implicit none
!Local variables
 integer  :: i,j,k
 real(8)  :: Q11,Q12,Q22,Q66,lower
 real(8),dimension(3,3) :: Ek_ply, Ek_xy,ply_trans
 real(8),dimension(NumLayers+1)::zk

! Initialization 
 Ae_ij = 0.d0
 Be_ij = 0.d0
 De_ij = 0.d0
 Q11 = 0.d0
 Q12 = 0.d0
 Q22 = 0.d0
 Q66 = 0.d0
 Ek_ply =0.d0
 Ek_xy = 0.d0
 ply_trans = 0.d0

!FlagIso = 0 => the layers are isotropic, else, orthotropic
 if (FlagIso == 0) then
   Q11 = 1.0
   Q22 = Q11
   Q12 = NU12
   Q66 = 0.5*(1-NU12)
 elseif (FlagIso == 1) then
   Q11 = E1/(1-NU12*NU21)
   Q22 = E2/(1-NU12*NU21)
   Q12 = (NU12*E2)/(1-NU12*NU21)
   Q66 = G12
 else
   print *, 'Error: FlagIso in nlams_abbd'
   stop
 endif

! In the case of the classical laminate theory, for the kth layer we obtain Eq. (14.17)
! The stress tensor at the point M is thus of the form shown in Eq. (14.18)
! The stress field reduces to the in-plane stresses, sigma_xx, sigma_yy, sigma_xy
! Ek_ply is composed of reudced stiffness matrix of kth layer 
! Refererence is book titled Composite matrials by Jean-Marie Berthelot
   Ek_ply(1,1) = Q11; Ek_ply(1,2) = Q12; Ek_ply(1,3) = 0.0
   Ek_ply(2,1) = Q12; Ek_ply(2,2) = Q22; Ek_ply(2,3) = 0.0
   Ek_ply(3,1) = 0.0; Ek_ply(3,2) = 0.0; Ek_ply(3,3) = Q66

! if(NumLayers == 0) then
! Ae_ij is same as the elasticity matrix (E) in local system
 
   Ae_ij = Ek_ply * E1/(1.d0-NU12*NU12)

!// Theoritically 
!  Ae_ij = Ek_ply * E1 * thickness /(1.d0-NU12*NU12)

! The membrane-bending coupling matrix B is zero when the laminate is symmetric,
! the symmetry implies a geometrical symmetry (thickness and z
! coordnates of their center) of their mechanical properies and
! of their orientations Berthelot's book entiled Composite Materials, p295
   Be_ij =  0.0d0 

! De_ij is bending stiffness matrix 
   De_ij = Ek_ply * E1 *(Thickness)**3/(12.d0*(1.d0-NU12*NU12))

! else

!// Composite model is available in different version of code

 end subroutine nlams_abbd

end module nlams_comp_abbd

