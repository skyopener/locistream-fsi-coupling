!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_TRANSFORMS. Satish Chimakurthi. 22Aug2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!     Contains several subroutines that compute transformation matrices
!
!-> Reference:
!->   Will be given seperately for each subroutine where applicable
!
!->Functions:
!
!     |-shell_trans
!
!
!->Remarks:
!    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_transforms
contains
  
  subroutine shell_transf(Xi,Xj,Xk,Yi,Yj,Yk,Zi,Zj,Zk,shell_trans)
    
    ! Notes:
    ! Computes transformation matrix from local frame to global frame
    ! For a triangle element with nodes i,j,k the local shell element frame is defined by fixing the local x-axis 
    ! along i-j with origin at node i, the local z-axis is defined as being normal to the plane of the element, and
    ! the local y-axis as the cross-product of the local x and z axes.
    
    ! Description of variables
    ! Xg - Nodal coordinates of an element in the global frame in the element's undeformed state
    
    ! I/O variables.
    real(8),intent(in)::Xi,Xj,Xk,Yi,Yj,Yk,Zi,Zj,Zk 
    real(8),dimension(3,3),intent(out)::shell_trans
    
    ! Local variables.
    real(8)::lji,lz
    
    Xji = Xj - Xi
    Yji = Yj - Yi
    Zji = Zj - Zi

    Xki = Xk - Xi
    Yki = Yk - Yi
    Zki = Zk - Zi

    XZijk = Yji*Zki - Zji*Yki
    YXijk = Zji*Xki - Xji*Zki
    ZYijk = Xji*Yki - Yji*Xki

    lji = sqrt((Xji**2) + (Yji**2) + (Zji**2))
    lz  = sqrt((XZijk)**2 + (YXijk**2) + (ZYijk**2))

    shell_trans(1,1) = Xji/lji
    shell_trans(1,2) = Yji/lji
    shell_trans(1,3) = Zji/lji

    shell_trans(3,1) = (1/lz)*(XZijk)
    shell_trans(3,2) = (1/lz)*(YXijk)
    shell_trans(3,3) = (1/lz)*(ZYijk)    
    
    shell_trans(2,1) = shell_trans(3,2)*shell_trans(1,3) - shell_trans(1,2)*shell_trans(3,3)
    shell_trans(2,2) = -shell_trans(3,1)*shell_trans(1,3) + shell_trans(1,1)*shell_trans(3,3)
    shell_trans(2,3) = shell_trans(3,1)*shell_trans(1,2) - shell_trans(1,1)*shell_trans(3,2)

!shell_trans at this point will transform the components of a vector in global frame to local frame

    shell_trans = transpose(shell_trans)

!shell_trans at this point will transform the components of a vector in local frame to global frame

  end subroutine shell_transf
  
  subroutine ply_transf(plyangle,ply_trans)
    
    ! Notes:
    ! Computes transformation matrix from global (fixed) frame to the local shell element frame.
    ! For a triangle element with nodes i,j,k the local shell element frame is defined by fixing the local x-axis
    ! along i-j with origin at node i, the local z-axis is defined as being normal to the plane of the element, and
    ! the local y-axis as the cross-product of the local x and z axes.
    
    ! Description of variables
    ! Xg - Nodal coordinates of an element in the global frame in the element's undeformed state
    
    ! I/O variables.
    real(8),intent(in)::plyangle
    real(8),dimension(3,3),intent(out)::ply_trans
    
    ! Local variables.
    real(8)::lji,lz
    
    ply_trans(1,1) = cos(plyangle)*cos(plyangle)
    ply_trans(1,2) = sin(plyangle)*sin(plyangle)
    ply_trans(1,3) = cos(plyangle)*sin(plyangle)
    
    ply_trans(2,1) = sin(plyangle)*sin(plyangle)
    ply_trans(2,2) = cos(plyangle)*cos(plyangle)
    ply_trans(2,3) = -cos(plyangle)*sin(plyangle)
    
    ply_trans(3,1) = -2*cos(plyangle)*sin(plyangle)
    ply_trans(3,2) = 2*cos(plyangle)*sin(plyangle)
    ply_trans(3,3) = cos(plyangle)*cos(plyangle) - sin(plyangle)*sin(plyangle)
    
  end subroutine ply_transf
  
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
end module nlams_transforms











