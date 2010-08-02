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
!     |-ply_trans
!     |-rodr_points2mat
!     |-vec_skew
!
!->Remarks:
! Computes transformation matrix from local frame to global frame: T_GE0
! For a triangle element with nodes i,j,k the local shell element frame is defined by fixing the local x-axis
! along i-j with origin at node i, the local z-axis is defined as being normal to the plane of the element, and
! the local y-axis as the cross-product of the local x and z axes.
!    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_transforms
! No global valiables

contains
  subroutine shell_transf(Xi,Xj,Xk,Yi,Yj,Yk,Zi,Zj,Zk,shell_trans)
    implicit none

! intent(in) 
    real(8),intent(in)::Xi,Xj,Xk,Yi,Yj,Yk,Zi,Zj,Zk 
! intent(out)
    real(8),dimension(3,3),intent(out)::shell_trans
    
! Local variables.
    real(8)::lji,lz
    real(8)::Xji,Yji,Zji,Xki,Yki,Zki,XZijk,YXijk,ZYijk
 !   real(8),dimension(3,3)::Temp,shell_trans_EG
 !   integer:: i

! Initialization of local variables
!    Xji  =0.0d0;Yji  =0.0d0;Zji  =0.0d0
!    Xki  =0.0d0;Yki  =0.0d0;Zki  =0.0d0
!    XZijk=0.0d0;YXijk=0.0d0;ZYijk=0.0d0
!    lji  =0.0d0;lz   =0.0d0  
!   
    Xji = Xj - Xi
    Yji = Yj - Yi
    Zji = Zj - Zi

    Xki = Xk - Xi
    Yki = Yk - Yi
    Zki = Zk - Zi

! Outer product based on (Xij,Yji,Zji) and (Xki,Yki,Zki)
    XZijk = Yji*Zki - Zji*Yki
    YXijk = Zji*Xki - Xji*Zki
    ZYijk = Xji*Yki - Yji*Xki

! Length of vectors (Xji,Yji,Zji) and (XZijk,YXijk,ZYijk)
    lji = sqrt((Xji*Xji) + (Yji*Yji) + (Zji*Zji))
    lz  = sqrt((XZijk*XZijk) + (YXijk*YXijk) + (ZYijk*ZYijk))

! unit vector along i-j axis
    shell_trans(1,1) = Xji/lji
    shell_trans(1,2) = Yji/lji
    shell_trans(1,3) = Zji/lji

! unit vecto normal to the undeformed element plane
    shell_trans(3,1) = XZijk/lz
    shell_trans(3,2) = YXijk/lz
    shell_trans(3,3) = ZYijk/lz    
    
    shell_trans(2,1) =  shell_trans(3,2)*shell_trans(1,3) - shell_trans(1,2)*shell_trans(3,3)
    shell_trans(2,2) = -shell_trans(3,1)*shell_trans(1,3) + shell_trans(1,1)*shell_trans(3,3)
    shell_trans(2,3) =  shell_trans(3,1)*shell_trans(1,2) - shell_trans(1,1)*shell_trans(3,2)

!    shell_trans_EG = shell_trans
! shell_trans at this point will transform the components of a vector in global frame to local frame: Transpose(T_GE0)

    shell_trans = transpose(shell_trans)

! shell_trans at this point will transform the components of a vector in local frame to global: T_GE0
    
!    Temp = matmul(shell_trans,shell_trans_EG)
  
!    do i=1,3
!      print *, Temp(i,1:3)
!    enddo
 end subroutine shell_transf
  
 subroutine ply_transf(plyangle,ply_trans)
    implicit none    
! I/O variables.
    real(8),intent(in)::plyangle
    real(8),dimension(3,3),intent(out)::ply_trans
    
! Local variables.
    real(8)::lji,lz
    
    ply_trans(1,1) = cos(plyangle)*cos(plyangle)
    ply_trans(1,2) = sin(plyangle)*sin(plyangle)
    ply_trans(1,3) = cos(plyangle)*sin(plyangle)
    
    ply_trans(2,1) =  sin(plyangle)*sin(plyangle)
    ply_trans(2,2) =  cos(plyangle)*cos(plyangle)
    ply_trans(2,3) = -cos(plyangle)*sin(plyangle)
    
    ply_trans(3,1) = -2.d0*cos(plyangle)*sin(plyangle)
    ply_trans(3,2) =  2.d0*cos(plyangle)*sin(plyangle)
    ply_trans(3,3) =  cos(plyangle)*cos(plyangle) - sin(plyangle)*sin(plyangle)
    
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











