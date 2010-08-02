!->Copyright by The University of Michigan, Aerospace Department. 2003
!
!->Module LIB_RODRIGUES. Rafa Palacios. 24Mar2003
!
!->Description.-
!
!  This module defines subroutines to operate with Rodrigues Parameters
!  (RodrPar) and rotation matrices to define transformations between
!  rectangular coordinates systems.
!
!->Subroutines:
!
!   rodr_add:         Compose two rotations defined by RodrPars.
!
!   rodr_beta2mat:    Convert cosines with b1 to Rotation Matrix.
!   rodr_dCdtheta:    Differential of a Rotation matrix with Theta.
!   rodr_dCTdtheta:   Differential of tranpose(Rot matrix) with Theta.
!   rodr_dCdtime:     Differential of a Rotation matrix with time.
!   rodr_dCTdtime:    Differential of tranpose(Rot matrix) with time.
!   rodr_dCdtdtheta:  Differential of a Rotation matrix with Time and Theta.
!   rodr_dCTdtdtheta: Differential of tranpose(Rot matrix) with Time and Theta.
!   rodr_dinvRdTheta: Differential of invR(Theta) with Theta.
!   rodr_dRdTheta:    Differential of R(Theta) with Theta.
!   rodr_dThetaAdd:   Differential of Theta02 w/r Theta01 or Theta12.
!
!   rodr_invrot:      Compute the inverse of the Rotational Operator R.
!   rodr_phi2mat:     Convert rotation angles to rotation matrix.
!   rodr_phi2theta:   Convert rotation angles to Rodrigues parameters.
!   rodr_points2mat:  Convert two points to Rotation Matrix.
!   rodr_mat2theta:   Convert a Rotation Matrix to RodrPars.
!   rodr_rotoperat:   Compute the Rotational Operator R(Theta).
!   rodr_theta2phi:   Convert RodrPars to simple rotation angles.
!   rodr_theta2mat:   Convert RodrPars to the Rotation Matrix.
!   rodr_rvec2mat:    Convert Cartesian rotation vector to rotation matrix.
!   rodr_rvec2rot:    Compute Rotational operator for a Cartesian rotation vector.
!
!   vect_cross:       Cross Product of two vectors.
!   vect_skew:        Convert a vector to its dual skew-symmetric matrix.
!
!-> Remarks.
!
!   1) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_rodrigues
 implicit none
!
!   (There are no public variables in this module).
!
!-> Private Variables.
!
  real(8), private, parameter,dimension(3,3) :: Unit= &
&       reshape((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
!
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_ADD.
!
!->Description.-
!
!   Compose two rotations defined by RodrPars.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_add (Theta01, Theta12)
!
! I/O Variables.
!
  real(8),intent(in)   :: Theta01 (3)
  real(8),intent(in)   :: Theta12 (3)   
  real(8),dimension(3) :: rodr_add ! Composed rotation (Theta02)
!
! Local variables.
!
  rodr_add= ( Theta01 +Theta12 -.5d0*vect_cross(Theta01,Theta12) ) &
&              / ( 1.d0 - .25d0 *dot_product(Theta01,Theta12))
!
  return
 end function rodr_add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subrutine RODR_BETA2MAT
! 
!-> Description.- 
! 
!  Find the transformation matrix from the director cosines for obliqueness.
!
!-> Remarks:
!
!  1) Definition of obliqueness follows the formulation of Popescu, Hodges
!     and Cesnik (2000): Beta is the vector of director cosines between the
!     the oblique axes and the beam longitudinal axis (normal axis to the
!     normal cross section, b1). They satisfy:
!
!               b1= Beta1*c1+Beta2*c2+Beta3*c3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rodr_beta2mat (Beta,RotMatrix)
!
!-> I/O Variables.
!
  real(8),intent(in) :: Beta(:)        ! Obliqueness Angles.
  real(8),intent(out):: RotMatrix(:,:) ! Rot matrix from b to c.
!
!-> Local variables.
!
  real(8):: Aux      ! Auxiliar Variable.
!
  Aux= 1.d0/(1.d0+Beta(1))
!
  RotMatrix(1,1)= Beta(1)
  RotMatrix(1,2)=-Beta(2)
  RotMatrix(1,3)=-Beta(3)
!
  RotMatrix(2,1)= Beta(2)
  RotMatrix(2,2)= 1.d0 - Beta(2)*Beta(2)*Aux
  RotMatrix(2,3)=-Beta(2)*Beta(3)*Aux
!
  RotMatrix(3,1)= Beta(3)
  RotMatrix(3,2)= RotMatrix(2,3)  
  RotMatrix(3,3)= 1.d0 - Beta(3)*Beta(3)*Aux
!
  return
 end subroutine rodr_beta2mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dCdTHETA.
!
!->Description.-
!
!   Differentiate a Rotation matrix C with respect to the 
!   components of the associated the Rodrigues Parameter.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dCdTheta (Theta)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)       ! RodrPars
  real(8), dimension(3,3,3) :: rodr_dCdTheta  ! Differential tensor.
!
!-> Local Variables.
!
  integer :: k              ! Counter
  real(8) :: Lambda         ! Defined as 1/(1+.25*Theta^2)
  real(8) :: SkewTheta(3,3) ! Dual Skew Symmetric form of Theta.
!
!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  SkewTheta=vect_skew(Theta)
!
! Compute the derivatives, component by component.
!
  do k=1,3
    rodr_dCdTheta(:,:,k)=                                                  &
&            Lambda*matmul(SkewTheta-Unit,vect_skew(Unit(:,k)))            &
&           +.5d0*Lambda*Lambda*(SkewTheta                                 &
&                               -.5d0*matmul(SkewTheta,SkewTheta))*Theta(k)
  end do
!
  return
 end function rodr_dCdTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dCTdTHETA.
!
!->Description.-
!
!   Differentiate a tranpose of a Rotation matrix C with respect to the
!   components of the associated the Rodrigues Parameter.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dCTdTheta (Theta)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)       ! RodrPars
  real(8), dimension(3,3,3) :: rodr_dCTdTheta  ! Differential tensor.
!
!-> Local Variables.
!
  integer :: k              ! Counter
  real(8) :: Lambda         ! Defined as 1/(1+.25*Theta^2)
  real(8) :: SkewTheta(3,3) ! Dual Skew Symmetric form of Theta.
!
!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  SkewTheta=vect_skew(Theta)
!
! Compute the derivatives, component by component.
!
  do k=1,3
    rodr_dCTdTheta(:,:,k)=                                                 &
&            Lambda*matmul(SkewTheta+Unit,vect_skew(Unit(:,k)))            &
&           -.5d0*Lambda*Lambda*(SkewTheta                                 &
&                               +.5d0*matmul(SkewTheta,SkewTheta))*Theta(k)
  end do
!
  return
 end function rodr_dCTdTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dCdTIME.
!
!->Description.-
!
!   Differentiate a Rotation matrix C with respect to the time
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dCdTime (Theta, Theta_t)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)       ! RodrPars
  real(8), intent(in)       :: Theta_t(3)     ! Time derivative of Theta.
  real(8), dimension(3,3)   :: rodr_dCdTime   ! Differentiated matrix.
!
!-> Local Variables.
!
  real(8) :: Lambda         ! Defined as 1/(1+.25*Theta^2)
  real(8) :: SkewTheta(3,3) ! Dual Skew Symmetric form of Theta.
!
!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  SkewTheta=vect_skew(Theta)
!
! Compute the derivatives.
!
    rodr_dCdTime= Lambda*(matmul(SkewTheta-Unit,vect_skew(Theta_t))) &
&  +.5d0*Lambda*Lambda*(SkewTheta - .5d0*matmul(SkewTheta,SkewTheta)) &
&                      *dot_product(Theta,Theta_t)
!
  return
 end function rodr_dCdTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dCTdTIME.
!
!->Description.-
!
!   Differentiate in time the transpose of a Rotation matrix, C.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dCTdTime (Theta,Theta_t)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)       ! RodrPars
  real(8), intent(in)       :: Theta_t(3)     ! Time derivative of Theta.
  real(8), dimension(3,3)   :: rodr_dCTdTime  ! Differentiated matrix.
!
!-> Local Variables.
!
  real(8) :: Lambda         ! Defined as 1/(1+.25*Theta^2)
  real(8) :: SkewTheta(3,3) ! Dual Skew Symmetric form of Theta.
!
!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  SkewTheta=vect_skew(Theta)
!
! Compute the derivatives.
!
    rodr_dCTdTime= Lambda*(matmul(SkewTheta+Unit,vect_skew(Theta_t))) &
&  -.5d0*Lambda*Lambda*(SkewTheta+.5d0*matmul(SkewTheta,SkewTheta))   &
&                      *dot_product(Theta,Theta_t)
!
  return
 end function rodr_dCTdTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dCdtdTheta
!
!->Description.-
!
!   Differentiate a rotation matrix, C, with respect to time and the cartesian 
!   components of its represenation using Rodrigues parameters, Theta.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dCdtdTheta (Theta,Theta_t)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)       ! RodrPars
  real(8), intent(in)       :: Theta_t(3)     ! Time derivative of Theta.
  real(8), dimension(3,3,3) :: rodr_dCdTdTheta! Differentiated matrix.
!
!-> Local Variables.
!
  integer :: k                ! Counter on the components of Theta.
  real(8) :: Lambda           ! Defined as 1/(1+.25*Theta^2)
  real(8) :: SkewTheta(3,3)   ! Dual Skew Symmetric form of Theta.
  real(8) :: SkewTheta_t(3,3) ! Dual Skew Symmetric form of dThetadT.
!
!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  SkewTheta=  vect_skew(Theta)
  SkewTheta_t=vect_skew(Theta_t)
!
! Compute the derivatives.
!
  do k=1,3
    rodr_dCdTdTheta(:,:,k)=                                              &
&       Lambda*matmul(SkewTheta_t,vect_skew(Unit(:,k)))                  &
&     +.5d0*Lambda*Lambda*(-matmul(SkewTheta-Unit,vect_skew(Unit(:,k)))  &
&                                 *dot_product(Theta,Theta_t)            &
&                          +(SkewTheta-.5d0*matmul(SkewTheta,SkewTheta)) &
&                          *(-Lambda*dot_product(Theta,Theta_t)*Theta(k) &
&                            +Theta_t(k))                                &
&                          +(SkewTheta_t-matmul(SkewTheta,SkewTheta_t))  &
&                          * Theta(k) )

  end do
!
  return
 end function rodr_dCdtdTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dCTdtdTheta
!
!->Description.-
!
!   Differentiate the transpose of a rotation matrix, C, with respect to time 
!   and the cartesian components of its represenation using Rodrigues 
!   parameters, Theta.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dCTdtdTheta (Theta,Theta_t)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)        ! RodrPars
  real(8), intent(in)       :: Theta_t(3)      ! Time derivative of Theta.
  real(8), dimension(3,3,3) :: rodr_dCTdtdTheta! Differentiated matrix.
!
!-> Local Variables.
!
  integer :: k                ! Counter on the components of Theta.
  real(8) :: Lambda           ! Defined as 1/(1+.25*Theta^2)
  real(8) :: SkewTheta(3,3)   ! Dual Skew Symmetric form of Theta.
  real(8) :: SkewTheta_t(3,3) ! Dual Skew Symmetric form of dThetadT.
!
!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  SkewTheta=  vect_skew(Theta)
  SkewTheta_t=vect_skew(Theta_t)
!
! Compute the derivatives.
!
  do k=1,3
    rodr_dCTdtdTheta(:,:,k)=                                             &
&       Lambda*matmul(SkewTheta_t,vect_skew(Unit(:,k)))                  &
&     -.5d0*Lambda*Lambda*( matmul(SkewTheta+Unit,vect_skew(Unit(:,k)))  &
&                         * dot_product(Theta,Theta_t)                   &
&                         + (SkewTheta+.5d0*matmul(SkewTheta,SkewTheta)) &
&                         * (Lambda*dot_product(Theta,Theta_t)*Theta(k)  &
&                            +Theta_t(k))                                &
&                         + (SkewTheta_t+matmul(SkewTheta,SkewTheta_t))  &
&                         *  Theta(k) )

  end do
!
  return
 end function rodr_dCTdtdTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dINVRdTHETA.
!
!->Description.-
!
!   Differentiate the inverse of the Rotation Operator (invR(Theta))
!   with respect to Theta.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dinvRdTheta (Theta)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)          ! RodrPars
  real(8), dimension(3,3,3) :: rodr_dinvRdTheta  ! Differential tensor.
!
!-> Local Variables.
!
  integer :: i,j,k         ! Counters.
!
! Compute the derivatives, component by component.
!
  do k=1,3
    rodr_dinvRdTheta(:,:,k)=.5d0*vect_skew(Unit(:,k))
    do i=1,3
      do j=1,3
        rodr_dinvRdTheta(i,j,k)= rodr_dinvRdTheta(i,j,k) +  &
&            .25d0*(Unit(i,k)*Theta(j)+Theta(i)*Unit(j,k))
      end do
    end do
  end do
!
  return
 end function rodr_dinvRdTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dRdTHETA.
!
!->Description.-
!
!   Differentiate the the Rotation Operator (R(Theta)) with Theta.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dRdTheta (Theta)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta(3)          ! RodrPars
  real(8), dimension(3,3,3) :: rodr_dRdTheta  ! Differential tensor.
!
!-> Local Variables.
!
  integer :: k              ! Counters.
  real(8) :: Lambda         ! Defined as 1/(1+.25*Theta^2)
  real(8) :: Unit(3,3) 
 
  Unit(1,1) = 1.0; Unit(2,2)=1.0; Unit(3,3) = 1.0

!
! Compute denominator (Lambda)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
!
! Compute the derivatives, component by component.
!
  do k=1,3
    rodr_dRdTheta(:,:,k)=-.5d0*Lambda*vect_skew(Unit(:,k))    &
&     -.5d0*Lambda*Lambda*(Unit-.5d0*vect_skew(Theta))*Theta(k)
  end do
!
  return
 end function rodr_dRdTheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_dThetaAdd
!
!->Description.-
!
!   Differentiate the compound rotation vector 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_dThetaAdd (Theta01,Theta12,Code)
!
!-> I/O Variables.
!
  real(8), intent(in)       :: Theta01(3)       ! RodrPars.
  real(8), intent(in)       :: Theta12(3)       ! Time derivative of Theta.
  integer, intent(in)       :: Code             ! =1: Derivative with Theta01.
                                                ! =2: Derivative with Theta12.    
  real(8), dimension(3,3)   :: rodr_dThetaAdd   ! Differentiated matrix.
!
!-> Local Variables.
!
  real(8) :: Aux33(3,3)        ! Auxiliary 3x3 matrix.
  real(8) :: BaseVec(3)        ! Unit Base Vector.
  integer :: i,j
  real(8) :: Theta02(3)        ! RodrPar vector of the compound rotation.
!
! Get the compound rotation.
!
  Theta02= rodr_add(Theta01,Theta12)
!
! Construct the numerator of the coefficients of the Jacobian matrix.
!
  Aux33=Unit

  do i=1,3
    BaseVec=0.d0
    BaseVec(i)=1.d0
    if (Code.eq.1) Aux33(i,:)= Aux33(i,:)-.5d0*vect_cross(BaseVec,Theta12)
    if (Code.eq.2) Aux33(i,:)= Aux33(i,:)-.5d0*vect_cross(Theta01,BaseVec)
!
    do j=1,3
      if (Code.eq.1) Aux33(i,j)=Aux33(i,j)+.25d0*Theta02(i)*Theta12(j)
      if (Code.eq.2) Aux33(i,j)=Aux33(i,j)+.25d0*Theta02(i)*Theta01(j)
    end do
  end do
!
! Add denominator.
!
  rodr_dThetaAdd= Aux33 / (1.d0-0.25d0*dot_product(Theta01,Theta12))
!
  return
 end function rodr_dThetaAdd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function RODR_PHI2MAT.
!
!-> Description.-
!
!  Convert simple rotations along the coordinate axis to rotation matrices.
!
!-> Remarks.
!
!  1) The elemental rotations are applied in the following order: Z,Y,X.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_phi2mat (Phi)
!
! I/O Variables.
!
  real(8),intent(in):: Phi(3)     ! Elemental rotations.
  real(8)::   rodr_phi2mat(3,3)   ! Rotation matrix.
!
! Local variables
!
  real(8):: RotX(3,3)      ! Simple rotations.
  real(8):: RotY(3,3)
  real(8):: RotZ(3,3)

  RotX(1,1)= 1.d0
  RotX(2,2)= cos(Phi(1))
  RotX(3,3)= cos(Phi(1))
  RotX(2,3)=-sin(Phi(1))
  RotX(3,2)= sin(Phi(1))
!
  RotY(2,2)= 1.d0
  RotY(3,3)= cos(Phi(2))
  RotY(1,1)= cos(Phi(2))
  RotY(3,1)=-sin(Phi(2))
  RotY(1,3)= sin(Phi(2))
!
  RotZ(3,3)= 1.d0
  RotZ(1,1)= cos(Phi(3))
  RotZ(2,2)= cos(Phi(3))
  RotZ(1,2)=-sin(Phi(3))
  RotZ(2,1)= sin(Phi(3))
!
  rodr_phi2mat=matmul(matmul(RotX,RotY),RotZ)
!
  return
 end function rodr_phi2mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function RODR_PHI2THETA.
!
!-> Description.-
!
!  Convert simple rotations along the coordinate axis to Rodrigues parameters.
!
!-> Remarks.
!
!  1) The elemental rotations are applied in the following order: Z,Y,X.
!
!  2) Input rotation (Phi) is given from the reference system (0) to the
!     rotated system (1)
!
!  3) The Rodrigues parameters are then Theta_0^1.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_phi2theta (Phi)
!
  real(8),intent(in):: Phi(3)    ! Elemental rotations.
  real(8):: rodr_phi2theta(3)    ! Rodrigues parameters.
  real(8):: RotMatrix(3,3)       ! Rotation matrix.
!
! Get the rotation matrix from 1 to 0 (C_01).
!
  RotMatrix=rodr_phi2mat(Phi)
!
! Get the Rodrigues parameters from 1 to 0 (Theta_0^1). 
!
  rodr_phi2theta=-rodr_mat2theta(RotMatrix)
!
  return
 end function rodr_phi2theta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function RODR_THETA2PHI.
!
!-> Description.-
!
!   Convert a RodrPars to Simple Rotations (rads) along the coordinate axis 
!   of the reference triad.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_theta2phi (Theta)
!
! I/O Variables.
!
  real(8),intent(in):: Theta(3)    ! RodrPars.
  real(8)::   rodr_theta2phi(3)    ! Rotation angles on each coordinate axis.
!
! Local variables
!
  integer :: i               ! Counter.
  real(8) :: ModTheta        ! Module of Theta
!
  ModTheta=sqrt(dot_product(Theta,Theta))
!
  if (ModTheta.eq.0.) then
    rodr_theta2phi=0.
  else
    do i=1,3
      rodr_theta2phi(i)=2.d0*atan(ModTheta/2.d0)*(Theta(i)/ModTheta)
    end do
  end if
!
  return
 end function rodr_theta2phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_THETA2MAT.
!
!->Description.-
!
!   Convert a RodrPars to the Rotation Matrix.
!
!-> Remarks.-
!
!   1) The rodrigues parameters from 1 to 2 (theta_1^2) correspond to the
!      rotation matrix from 1 to 2: C^2^1.
!
!      A vector given in its components in the base 1 (u_1), is expressed in its
!      components in a base 2 (u_2) as
!
!                  u_2= C^2^1*u_1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_theta2mat (Theta)
!
! I/O Variables.
!
  real(8),intent(in)    :: Theta(3)        ! RodrPars.
  real(8),dimension(3,3):: rodr_theta2mat  ! Rotation Matrix.
!
! Local variables
!
  real(8) :: RotMatrix(3,3) ! Rotation Matrix.
  real(8) :: Denom          ! Denominator.
!
  Denom=1+(Theta(1)**2+Theta(2)**2+Theta(3)**2)*0.25d0
!
  RotMatrix(1,1)=(2-Denom)+0.5*Theta(1)**2
  RotMatrix(2,2)=(2-Denom)+0.5*Theta(2)**2
  RotMatrix(3,3)=(2-Denom)+0.5*Theta(3)**2
!
  RotMatrix(1,2)= Theta(3)+0.5*Theta(2)*Theta(1)
  RotMatrix(1,3)=-Theta(2)+0.5*Theta(3)*Theta(1)
  RotMatrix(2,1)=-Theta(3)+0.5*Theta(2)*Theta(1)
  RotMatrix(2,3)= Theta(1)+0.5*Theta(2)*Theta(3)
  RotMatrix(3,1)= Theta(2)+0.5*Theta(3)*Theta(1)
  RotMatrix(3,2)=-Theta(1)+0.5*Theta(2)*Theta(3)
!
  rodr_theta2mat=RotMatrix/Denom
!
  return
 end function rodr_theta2mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_RVEC2MAT.
!
!->Description.-
!
!   Convert Cartesian rotation vector to rotation matrix form.
!
!-> Remarks.-
!
!   1) The cartesian rotation vector from 1 to 2 (alpha_1^2) corresponds to the
!      rotation matrix from 1 to 2: C^2^1.
!
!      A vector given in its components in the base 1 (u_1), is expressed in its
!      components in a base 2 (u_2) as
!
!                  u_2= C^2^1*u_1
!
!  2) Equation (4.5) of (Cardona and Geradin, 2001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_rvec2mat (Alpha)

! I/O Variables.
  real(8),intent(in)    :: Alpha(3)       ! Cartesian rotation vector.
  real(8),dimension(3,3):: rodr_rvec2mat  ! Rotation Matrix.

! Local variables
  real(8):: AlphaHat(3,3)  ! Rotation Matrix.
  real(8):: Norm           ! Denominator.
  real(8):: RotMatrix(3,3) ! Rotation Matrix.
!
  Norm= sqrt(Alpha(1)**2+Alpha(2)**2+Alpha(3)**2)

  if (Norm.gt.0.d0) then
    AlphaHat=vect_skew(Alpha)
    RotMatrix= Unit &
&            + sin(Norm)/Norm * AlphaHat &
&            + (1.d0-cos(Norm))/(Norm*Norm) * matmul(AlphaHat,AlphaHat)
  else
    RotMatrix= Unit
  end if

  rodr_rvec2mat=transpose(RotMatrix)
  return
 end function rodr_rvec2mat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_RVEC2ROT.
!
!->Description.-
!
!   Obtain rotational (tangential) operator from Cartesian rotation vector.
!
!-> Remarks.-
!
!   1) The cartesian rotation vector from 1 to 2 (alpha_1^2) corresponds to the
!      rotation matrix from 1 to 2: C^2^1.
!
!      A vector given in its components in the base 1 (u_1), is expressed in its
!      components in a base 2 (u_2) as
!
!                  u_2= C^2^1*u_1
!
!  2) Equation (4.11) of (Cardona and Geradin, 2001).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_rvec2rot (Alpha)

! I/O Variables.
  real(8),intent(in)    :: Alpha(3)       ! Cartesian rotation vector.
  real(8),dimension(3,3):: rodr_rvec2rot  ! Rotational operator.

! Local variables
  real(8):: AlphaHat(3,3)  ! Rotation Matrix.
  real(8):: Norm           ! Denominator.
  real(8):: RotMatrix(3,3) ! Rotation Matrix.
!
  Norm= sqrt(Alpha(1)**2+Alpha(2)**2+Alpha(3)**2)

  if (Norm.gt.0.d0) then
    AlphaHat=vect_skew(Alpha)
    RotMatrix= Unit &
&            + (cos(Norm)-1.d0)/(Norm*Norm) * AlphaHat &
&            + ((1.d0-cos(Norm)/Norm)/(Norm*Norm)) * matmul(AlphaHat,AlphaHat)
  else
    RotMatrix= Unit
  end if

  rodr_rvec2rot=transpose(RotMatrix)
  return
 end function rodr_rvec2rot




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine RODR_POINTS2MAT.
!
!->Description.-
!
!   Convert the location of two points in frame A to a Rotation Matrix from A to B.
!
!   Point 1 is in the Z axis of B.
!   Point 2 is in the X-Z plane of B.
!
!->Remarks.-
!
!   1) Both coordinate reference systems have the same origin, 0.
!
!   2) The rotation matrix from A to B is defined as:
!
!          Bi=(C_BA)ij·Aj
!
!          Bi and Aj are the base vectors of frames B and A, respectively.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rodr_points2mat (Point1, Point2, RotMatrix)
!
!-> I/O Variables.
!
  real(8),intent(in) :: Point1(3)       ! Point 1 coordinates.
  real(8),intent(in) :: Point2(3)       ! Point 2 coordinates.
  real(8),intent(out):: RotMatrix(3,3)  ! Rotation Matrix.
!
!-> Local Variables.
!
  real(8) :: b1(3)    ! Unit vector of B along x1.
  real(8) :: b2(3)    ! Unit vector of B along x2.
  real(8) :: b3(3)    ! Unit vector of B along x3.
  real(8) :: Unit1(3) ! Unit vector in the direction of point 1.
  real(8) :: Unit2(3) ! Unit vector in the direction of point 2.
!
! Get unit vectors in the direction of point 1 and point2.
!
  Unit1=Point1/dsqrt(dot_product(Point1,Point1))
  Unit2=Point2/dsqrt(dot_product(Point2,Point2))
!
! Vector 0-1 gives b3.
!
  b3= Unit1
!
! b2 is normal to 0-1 and 0-2.
!
  b2= matmul(vect_skew(Unit1),Unit2)
  b2= b2/dot_product(b2,b2)
!
! b1 is normal to b2 and b3.
!
  b1= matmul(vect_skew(b2),b3)
!
! The Rotation matrix from A to B is defined by the components of this
! orthonormal base.
!
  RotMatrix(1,:)=b1
  RotMatrix(2,:)=b2
  RotMatrix(3,:)=b3
!  
  return
 end subroutine rodr_points2mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function RODR_MAT2THETA.

!-> Description.-

!   Convert a Rotation Matrix to RodrPars.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_mat2theta (RotMatrix)

! I/O Variables.

  real(8),intent(in)  :: RotMatrix(3,3)    ! Rotation Matrix
  real(8),dimension(3):: rodr_mat2theta(3) ! RodrPars.

! Local variables

  integer  :: i               ! Counter.
  real(8)  :: ThetaSkew(3,3)  ! Skew-Symmetric Matrix from RodrPars.
  real(8)  :: Trace           ! Trace of the RotMatrix.

  Trace=0.d0
  do i=1,3
    Trace=Trace+RotMatrix(i,i)
  end do

  ThetaSkew= -2.d0*(RotMatrix-transpose(RotMatrix))/(1.d0+Trace)

  rodr_mat2theta(1)= ThetaSkew(3,2)
  rodr_mat2theta(2)= ThetaSkew(1,3)
  rodr_mat2theta(3)= ThetaSkew(2,1)

  return
 end function rodr_mat2theta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function RODR_ROTOPERAT.
!
!->Description.-
!
!   Compute the Rotational Operator R(Theta).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_rotoperat (Theta)
!
! I/O Variables.
!
  real(8),intent(in) ::    Theta(3)        ! RodrPars
  real(8),dimension(3,3):: rodr_rotoperat  ! Rotational Operator.
!
! Local variables
!
  real(8) :: Lambda         ! Defined as 1/(1+.25*Theta^2)
!
  Lambda=1.d0/(1.d0+.25d0*dot_product(Theta,Theta))
  rodr_rotoperat=Lambda*(Unit-.5d0*vect_skew(Theta))
!
  return
 end function rodr_rotoperat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Function RODR_INVROT.
!
!->Description.-
!
!   Compute the Inverse of the Rotational Operator R(Theta).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function rodr_invrot (Theta)
!
! I/O Variables.
!
  real(8),intent(in)     :: Theta(3)    ! RodrPars
  real(8),dimension(3,3) :: rodr_invrot ! Inverse Rotational Operator.
!
! Local variables
!
  integer :: i,j             ! Counter.
  real(8) :: ThetaSkew(3,3)  ! Skew-Symmetric Matrix from RodrPars.
!
  ThetaSkew=vect_skew(Theta) 
!
  do i=1,3
    do j=1,3
      rodr_invrot(i,j)=Unit(i,j)+.5d0*ThetaSkew(i,j)+.25d0*Theta(i)*Theta(j)
    end do
  end do
!
  return
 end function rodr_invrot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine VECT_CROSS.
!
!->Description.-
!
!   Cross Product of two vectors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function vect_cross (Vec1, Vec2)
!
! I/O Variables.
!
  real(8),intent(in)  :: Vec1 (3)
  real(8),intent(in)  :: Vec2 (3)
  real(8),dimension(3):: vect_cross
!
! Local variables.
!
  vect_cross(1)= Vec1(2)*Vec2(3) - Vec1(3)*Vec2(2)
  vect_cross(2)= Vec1(3)*Vec2(1) - Vec1(1)*Vec2(3)
  vect_cross(3)= Vec1(1)*Vec2(2) - Vec1(2)*Vec2(1)
!
  return
 end function vect_cross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine VECT_SKEW.
!
!->Description.-
!
!   Compute the Skew-Symmetric Matrix given a vector.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function vect_skew (Vector)
!
!-> I/O Variables.
!
  real(8), intent(in) :: Vector(3)
  real(8), dimension(3,3):: vect_skew
!
  vect_skew=0.0d0
  vect_skew(1,2)=-Vector(3)
  vect_skew(1,3)= Vector(2)
  vect_skew(2,1)= Vector(3)
  vect_skew(2,3)=-Vector(1)
  vect_skew(3,1)=-Vector(2)
  vect_skew(3,2)= Vector(1)
!
  return
 end function vect_skew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_rodrigues
