module lib_fem
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine FEM_GAUSS_VAL
!
!-> Description.-
!
!   Database of Weights and Coordinates for Gaussian Quadrature.
!
!->Remarks:
!
!   1) Coefficients are extracted from Zienkiewitz' book "The Finite
!      Element Method".
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine fem_gauss_val (Flag, NumGauss, Coords, Weight)
!
 use lib_kind
!-> I/O Variables.
!
  integer,intent(in) :: Flag              ! =0 for CTRIA, /=0 for CQUAD.
  integer,intent(in) :: NumGauss          ! Number of Gaussian points.
  real(8),intent(out):: Coords(:,:)       ! Local Coordinates of the Points.
  real(8),intent(out):: Weight(:)         ! Weights of the Gaussian Points.
!
!-> Local variables.
!
  integer, parameter:: MaxPoints=10 ! Maximum # of linear points.
  integer :: iGauss                 ! Counter in the Gaussian points.
  integer :: i1,i2                  ! Counters in the linear points.
  integer :: NumPoints              ! Number of linear points.
  real(8) :: CoordsLin(MaxPoints)   ! Local Coordinates for [-1,1] integration.
  real(8) :: WeightLin(MaxPoints)   ! Local Coordinates for [-1,1] integration.
!
  
! Default output values.
!
  Coords=0.d0
  Weight=0.d0
!
! Quadrilateral Elements.
!
  if (Flag.ne.0) then
    NumPoints=nint(sqrt(float(NumGauss)))
    if ((NumGauss-NumPoints*NumPoints).ne.0) NumPoints=0
    select case(NumPoints)
    case(1)
      CoordsLin(1)=0.d0
      WeightLin(1)=0.d0
    case(2)
      CoordsLin(1)=-1.d0/dsqrt(3.d0)
      CoordsLin(2)= 1.d0/dsqrt(3.d0)
      WeightLin(1)= 1.d0
      WeightLin(2)= 1.d0
    case(3)
      CoordsLin(1)=-dsqrt(.6d0)
      CoordsLin(2)= 0.d0
      CoordsLin(3)= dsqrt(.6d0)
      WeightLin(1)= 5.d0/9.d0
      WeightLin(2)= 8.d0/9.d0
      WeightLin(3)= 5.d0/9.d0
    case(4)
      CoordsLin(1)=-0.861136311594953
      CoordsLin(2)=-0.339981043584856
      CoordsLin(3)= 0.339981043584856
      CoordsLin(4)= 0.861136311594953
      WeightLin(1)= 0.347854845137454
      WeightLin(2)= 0.652145154862546
      WeightLin(3)= 0.652145154862546
      WeightLin(4)= 0.347854845137454
    case(5)
      CoordsLin(1)=-0.906179845938664
      CoordsLin(2)=-0.538469310105683
      CoordsLin(3)= 0.d0
      CoordsLin(4)= 0.538469310105683
      CoordsLin(5)= 0.906179845938664
      WeightLin(1)= 0.236926885056189
      WeightLin(2)= 0.478628670499366
      WeightLin(3)= 0.568888888888889
      WeightLin(4)= 0.478628670499366
      WeightLin(5)= 0.236926885056189
    case default
    end select
!
    iGauss=0
    do i1=1,NumPoints
      do i2=1,NumPoints
        iGauss= iGauss+1
        Coords(iGauss,1)= CoordsLin(i1)
        Coords(iGauss,2)= CoordsLin(i2)
        Weight(iGauss)  = WeightLin(i1)*WeightLin(i2)
      end do
    end do
!
! Triangular elements.
!
  else if (Flag.eq.0) then
    NumPoints=NumGauss

    select case(NumPoints)
    case(3)
      Coords(1,1)=1.d0/6.d0
      Coords(1,2)=1.d0/6.d0
      Coords(2,1)=2.d0/3.d0
      Coords(2,2)=1.d0/6.d0
      Coords(3,1)=1.d0/6.d0
      Coords(3,2)=2.d0/3.d0
      Weight(1:3)=1.d0/6.d0
    case(7)
      Coords(1,1)= 1.0d0/3.0d0
      Coords(1,2)= 1.0d0/3.0d0
      Weight(1)  = 0.225000000000000
      Coords(2,1)= 0.470142064105115
      Coords(2,2)= 0.470142064105115
      Weight(2)  = 0.132394152788506
      Coords(3,1)= 0.059715871789770
      Coords(3,2)= 0.470142064105115
      Weight(3)  = 0.132394152788506
      Coords(4,1)= 0.470142064105115
      Coords(4,2)= 0.059715871789770
      Weight(4)  = 0.132394152788506
      Coords(5,1)= 0.101286507323456
      Coords(5,2)= 0.101286507323456
      Weight(5)  = 0.125939180544827
      Coords(6,1)= 0.797426985353087
      Coords(6,2)= 0.101286507323456
      Weight(6)  = 0.125939180544827
      Coords(7,1)= 0.101286507323456
      Coords(7,2)= 0.797426985353087
      Weight(7)  = 0.125939180544827
    case default
end select
  end if
!
  return
 end subroutine fem_gauss_val
end module lib_fem
