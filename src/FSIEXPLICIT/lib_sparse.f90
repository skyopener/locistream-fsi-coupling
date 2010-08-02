!->Copyright by The University of Michigan, Aerospace Department. 2006
!
!->Module LIB_SPARSE. Rafa Palacios. 30Aug2006
!
!->Description.-
!
!  This module defines tools for handling sparse matrix.
!
!->Subroutines:
!
!   sparse_addmat:    Add submatrix to sparse matrix.
!   sparse_addval:    Add element to sparse matrix.
!   sparse_bandwidth: Compute bandwidth of sparse matrix.
!   sparse_copy2band: Write sparse matrix in banded form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_sparse
 implicit none
!
! Define derived types.
!
  type sparse
    integer i
    integer j
    real(8) a
  end type

 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_ADDVAL.
!
!->Description.-
!
!   Add element to sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_addval (i,j,Val,DimArray,Matrix)

  integer,     intent(in)   :: i, j      ! Coefficient in the matrix.
  real(8),     intent(in)   :: Val       ! Value to be added.
  integer,     intent(inout):: DimArray  ! Current storage dimension.
  type(sparse),intent(inout):: Matrix(:) ! Sparse matrix.

  logical:: Flag   ! Logical flag.
  integer:: k      ! Counter of element

! If there is already a (i,j) term in the matrix, add the new value.

  Flag=.false.
  do k=1,DimArray
    if ((Matrix(k)%i.eq.i).and.(Matrix(k)%j.eq.j)) then
      Matrix(k)%a=Matrix(k)%a+Val
      Flag=.true.
      exit
    end if
  end do

! If (i,j) was empty, create a new entry.

  if (.not.Flag) then

    if (DimArray.eq.size(Matrix)) then
      stop 'ERROR: Not enough memory for sparse matrix allocation'
    end if

    DimArray=DimArray+1
    Matrix(DimArray)%i= i
    Matrix(DimArray)%j= j
    Matrix(DimArray)%a= Val
  end if

  return
 end subroutine sparse_addval


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_ADDMAT.
!
!->Description.-
!
!   Add submatrix to sparse matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_addmat (i1,j1,Mat,DimArray,SprMat)

  integer,     intent(in)   :: i1,j1     ! First coefficient on the SprMat.
  real(8),     intent(in)   :: Mat(:,:)  ! Value to be added.
  integer,     intent(inout):: DimArray  ! Current storage dimension.
  type(sparse),intent(inout):: SprMat(:) ! Sparse matrix.

  integer:: i,j    ! Counter on the element of the matrix.

  do i=1,size(Mat,DIM=1)
    do j=1,size(Mat,DIM=2)
      if (Mat(i,j).ne.0.d0) then
        call sparse_addval (i1+i,j1+j,Mat(i,j),DimArray,SprMat)
      end if
    end do
  end do

  return
 end subroutine sparse_addmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_GETBANDWIDTH.
!
!->Description.-
!
!   Get lower and upper bandwidths of sparse matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_getbandwidth (DimArray,Matrix,KL,KU)

  integer,     intent(in) :: DimArray  ! Current storage dimension.
  type(sparse),intent(in) :: Matrix(:) ! Sparse matrix.
  integer,     intent(out):: KL        ! Number of subdiagonals of the banded form.
  integer,     intent(out):: KU        ! Number of superdiagonals of the banded form.

  integer :: k  ! Counter on the elements of the array.

  KL=0
  KU=0

  do k=1,DimArray
    if ((Matrix(k)%i-Matrix(k)%j).gt.KL) KL=Matrix(k)%i-Matrix(k)%j
    if ((Matrix(k)%j-Matrix(k)%i).gt.KU) KU=Matrix(k)%j-Matrix(k)%i
  end do

  return
 end subroutine sparse_getbandwidth


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine SPARSE_COPY2BAND.
!
!->Description.-
!
!   Copy matrix to bandwidth form used in LAPACK.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine sparse_copy2band (DimArray,Matrix,KL,KU,BandedMat)

  integer,     intent(in) :: DimArray       ! Current storage dimension.
  type(sparse),intent(in) :: Matrix(:)      ! Sparse matrix.
  integer,     intent(in) :: KL             ! Number of superdiagonals of the banded form.
  integer,     intent(in) :: KU             ! Number of subdiagonals of the banded form.
  real(8),     intent(out):: BandedMat(:,:) ! Number of superdiagonals of the banded form.

  integer :: k  ! Counter on the elements of the array.

  do k=1,DimArray
    BandedMat(KL+KU+1+Matrix(k)%i-Matrix(k)%j,Matrix(k)%j) = Matrix(k)%a
  end do

  return
 end subroutine sparse_copy2band


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_sparse
