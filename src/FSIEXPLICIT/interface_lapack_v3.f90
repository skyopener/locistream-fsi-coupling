!->Copyright by The University of Michigan, Aerospace Department. 2006
!
!->Module INTERFACE_LAPACK. Rafa Palacios. 04Oct2006
!
!->Description.-
!
!  This module includes interface for real*8 LAPACK Version 3 routines.
!
!->Subroutines:
!
!     -lapack_inv:         Inverse of general matrix.
!     -lapack_lufact:      LU factorization of general matrix.
!     -lapack_luback:      LU backsubstution of general matrix.
!     -lapack_band_lufact: LU factorization of banded matrix.
!     -lapack_band_inv:    Inverse of banded matrix.
!     -lapack_nonsymeigv:  Complex eigenvalues of generalized eigv problem.
!
!->Remarks.-
!
!  1) External routines are part of the LAPACK 3 library, which needs
!     to be compiled with this module.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_lapack
 implicit none

!   (There are no public variables in this module).

 contains

! NEW SUBROUTINE
 subroutine inv_matrix(N,A,InvA,Info)
 implicit none
 integer,intent(in)::     N
 real(8),intent(in)::     A(:,:)
 real(8),intent(out):: InvA(:,:)
 integer,intent(out):: Info

! local variables
 real(8):: max,eps
 real(8):: w
 integer:: p,q
 integer:: i,j,k,iw,jw
 integer:: M

 integer,allocatable::iwork(:)
 integer,allocatable::jwork(:)
 real(8),allocatable::B(:,:)

 allocate(B(N,N*2),iwork(N),jwork(N))

 M=2*N
 eps=1.d-32

 do j=1,M
  do i=1,N
   if(j.le.N)then
    B(i,j)=A(i,j)
   else
    B(i,j)=0.0d0
    if(i.eq.j-N) B(i,j)=1.0d0
   endif
  enddo
 enddo

 do i=1,N
  iwork(i)=i
  jwork(i)=i
 enddo

! Do loop
 do k=1,N
  max = dabs(B(k,k))
  p   = k
  q   = k
   do j=k,N
   do i=k,N
    if(max.lt.dabs(B(i,j))) then
     max = dabs(B(i,j))
     p   = i
     q   = j
    endif
   enddo
   enddo
   if(max.gt.eps)then
    do i=1,N
     w     = B(i,k)
     B(i,k)= B(i,q)
     B(i,q)= w
    enddo
    do j=k,N+k-1
     w     = B(k,j)
     B(k,j)= B(p,j)
     B(p,j)= w
    enddo

    i       = iwork(k)
    iwork(k)= iwork(p)
    iwork(p)= i
    j       = jwork(k)
    jwork(k)= jwork(p)
    jwork(p)= j

    do j=k+1,M
     B(k,j) = B(k,j)/B(k,k)
    enddo

    do i=1,N
     if(i.ne.k) then
      do j=k+1,M
       B(i,j) = B(i,j) - B(i,k)*B(k,j)
      enddo
     endif
    enddo
   else ! max.gt.eps
    print *, '*** Matrix is ill ***'
    Info = 1
   endif
  enddo ! k-loop

  do j=N+1,M
   do i=1,N
    iw     = jwork(i)
    B(iw,N)= B(i,j)
   enddo
   do i=1,N
    B(i,j) = B(i,N)
   enddo
  enddo

  do i=1,N
   do j=N+1,M
    jw     = iwork(j-N)
    B(N,jw)= B(i,j)
   enddo
   do j=N+1,M
    B(i,j) = B(N,j-N)
   enddo
  enddo

! Solution
  do j=1,N
  do i=1,N
   InvA(i,j)=B(i,j+N)
  enddo
  enddo

  deallocate(B,iwork,jwork)

  return
  end subroutine inv_matrix

  subroutine lu_decomp(N, A, IP, Info)

  implicit none
  integer,intent(in)::     N
  real(8),intent(inout)::  A(:,:)
  integer,intent(out):: Info
  integer,intent(out):: IP(:)
  real(8),allocatable:: W(:)
!  integer,allocatable:: IP(:)

! local variables
  integer :: i,j,k
  integer :: L,LV
  real(8) :: eps
  real(8) :: AL

  allocate(W(N))

  eps = 1.0e-30

  do i=1,N
   IP(i) = i
  enddo

  do k=1,N
   L  = k
   AL = dabs(A(IP(L),k))
    do i=k+1,N
     if(dabs(A(IP(i),k)).gt.AL)then
      L  = i
      AL = dabs(A(IP(L),k))
     endif
    enddo

    if(L.ne.k)then
! --- Row exchange ---
     LV    = IP(k)
     IP(k) = IP(L)
     IP(L) = LV
    endif

    if(dabs(A(IP(k),k)).gt.eps) then
! --- Gauss elimination ---
     A(IP(k),k) = 1.0d0 / A(IP(k),k)

     do i= k+1, N
      A(IP(i),k) = A(IP(i),k) * A(IP(k),k)
      do j= k+1, N
       W(j) = A(IP(i),j) - A(IP(i),k) * A(IP(k),j)
      enddo
      do j= k+1, N
       A(IP(i),j) = W(j)
      enddo
     enddo
    else
     print *, '*** Matrix singular in lu_decomp ***'
     print *, '*** -th pivot ***',k
     info = 1
    endif
  enddo
   
!  do j=1,N
!  do i=1,N
!    B(i,j) = A(i,j)
!  enddo
!  enddo
 
  deallocate(W)
   
  return
  end subroutine lu_decomp
     
  subroutine lu_solve(N,A,IP,B,X)

  implicit none

  integer, intent(in) :: N
  integer, intent(in) :: IP(:)
  real(8), intent(in) :: B(:)
  real(8), intent(in) :: A(:,:)
  real(8), intent(inout):: X(:)
  
  integer :: i,j,k
  real(8) :: T
  
  T=0.d0

!-----Forward substitution ----
   do i=1, N
    T   = B(IP(i))
    do j=1, i-1
     T  = T - A(IP(i),j) * real(X(j))
    enddo
    X(i)= T
   enddo  

!-----Backward substitution ----
   do i=N, 1, -1
    T = X(i)
    do j=i+1, N
     T = T - A(IP(i),j) * real(X(j))
    enddo
    X(i) = T * A(IP(i),i)
   enddo

   return
   end subroutine lu_solve
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_INV
!
!->Description.-
!
!   Computes the inverse of a general matrix, using the LU factorization
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrf:
!        Computes an LU factorization of a general matrix, using partial
!        pivoting with row interchanges.
!
!      -dgetrI:
!         Computes the inverse of a general matrix, using the LU factorization
!         computed by DGETRF.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_inv (N,A,invA,Info)

  integer,intent(in) :: N         ! Size of the matrix.
  real(8),intent(in) :: A   (:,:) ! Matrix to be inverted.
  real(8),intent(out):: invA(:,:) ! A^-1.
  integer,intent(out):: Info      ! Error codes.

  integer,allocatable:: LUpivot(:)
  real(8),allocatable:: LUwork (:)

! Allocate auxiliary memory.

  allocate (LUpivot(N));  LUpivot=0
  allocate (LUwork(N*N)); LUwork =0.d0

! Perform LU decomposition.
!  do i=1,N
!  do j=1,N
!    invA(i,j)=A(i,j)
!  enddo
!  enddo
!  print *, 'Before DGETRF',N
  invA=A
  call DGETRF(N,N,invA(1:N,1:N),N,LUpivot(1:N),Info)
  if (Info.ne.0) return
!  print *, 'Info:',Info

! Perform matrix inversion.
   call DGETRI(N,invA(1:N,1:N),N,LUpivot(1:N),LUwork(1:N*N),N*N,Info)
!  print *, 'Info:',Info
  
  deallocate (LUpivot,LUwork)
  return
 end subroutine lapack_inv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_LUFACT
!
!->Description.-
!
!   Computes the LU factorization of a general matrix.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrf:
!        Computes an LU factorization of a general matrix, using partial
!        pivoting with row interchanges.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_lufact (N,A,LUpivot,Info)

  integer,intent(in)   :: N         ! Size of the matrix.
  real(8),intent(inout):: A   (:,:) ! LU factorization of A.
  integer,intent(out)  :: LUpivot(:)! LU pivoting index.
  integer,intent(out)  :: Info      ! Error codes.

! Perform LU decomposition.

  call DGETRF (N,N,A(1:N,1:N),N,LUpivot(1:N),Info)

  return
 end subroutine lapack_lufact


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_LUBACK
!
!->Description.-
!
!   Solves a general system of linear equations AX=B.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrs:
!        Solves a general system of linear equations AX=B, A**T X=B
!        or A**H X=B, using the LU factorization computed by DGETRF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_luback (N,A,LUpivot,B,Info)

  integer,intent(in)   :: N            ! Size of the matrix.
  real(8),intent(in)   :: A      (:,:) ! Matrix and its LU factorization.
  integer,intent(in)   :: LUpivot(:)   ! LU pivoting index.
  real(8),intent(inout):: B      (:)   ! B in input, X in output.
  integer,intent(out)  :: Info         ! Error codes.

  real(8),allocatable  :: X(:,:)

  allocate (X(N,1))
  X(1:N,1)=B(1:N)

! Solve AX=B.

  call DGETRS ('N',N,1,A(1:N,1:N),N,LUpivot(1:N),X(1:N,1:1),N,Info)

  B(1:N)=X(1:N,1)

  deallocate (X)
  return
 end subroutine lapack_luback


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_BAND_LUFACT
!
!->Description.-
!
!   Computes the LU factorization of a banded matrix.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgetrf:
!        Computes an LU factorization of a general matrix, using partial
!        pivoting with row interchanges.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_band_lufact (N,KL,KU,A,LUpivot,Info)

  integer,intent(in)   :: N         ! Size of the matrix.
  integer,intent(in)   :: KL        ! Number of subdiagonals in the band of B.
  integer,intent(in)   :: KU        ! Number of superdiagonals in the band of B.
  real(8),intent(inout):: A   (:,:) ! A; LU factorization of A.
  integer,intent(out)  :: LUpivot(:)! LU pivoting index.
  integer,intent(out)  :: Info      ! Error codes.

! Perform LU decomposition.

  call DGBTRF (N,N,KL,KU,A(1:2*KL+KU+1,1:N),2*KL+KU+1,LUpivot(1:N),Info)

  return
 end subroutine lapack_band_lufact


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_BAND_INV
!
!->Description.-
!
!   Computes the inverse of a banded matrix, using the LU factorization
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -dgbtrf:
!        Computes an LU factorization of a banded matrix, using partial
!        pivoting with row interchanges.
!
!      -dgbtri:
!         Computes the inverse of a banded matrix, using the LU factorization
!         computed by DGBTRF.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine lapack_band_inv (N,KL,KU,A,invA,Info)

  integer,intent(in) :: N         ! Size of the matrix.
  integer,intent(in) :: KL        ! Number of subdiagonals in the band of B.
  integer,intent(in) :: KU        ! Number of superdiagonals in the band of B.
  real(8),intent(in) :: A   (:,:) ! Banded atrix to be inverted.
  real(8),intent(out):: invA(:,:) ! A^-1 (full matrix)
  integer,intent(out):: Info      ! Error codes.

! Local variables.

  integer:: i
  real(8),allocatable:: LUdecom(:,:)
  integer,allocatable:: LUpivot(:)
  real(8),allocatable:: LUwork (:)

! Allocate auxiliary memory.

  allocate (LUpivot(N));  LUpivot=0
  allocate (LUwork(N*N)); LUwork =0.d0
  allocate (LUdecom(2*KL+KU+1,N)); LUdecom=A

! Perform LU decomposition.

  call DGBTRF (N,N,KL,KU,LUdecom,2*KL+KU+1,LUpivot,Info)
  if (Info.ne.0) return

! Perform matrix inversion.

  invA=0.d0
  do i=1,N
    invA(i,i)=1.d0
    call DGBTRS ('N',N,KL,KU,1,LUdecom,2*KL+KU+1,LUpivot,invA(i,:),N,Info)
  end do

  deallocate (LUpivot,LUwork,LUdecom)
  return
 end subroutine lapack_band_inv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine LAPACK_NONSYMEIGV
!
!->Description.-
!
!   Computes the complex eigenvalues of the generalized problem AX= lambda*B*X.
!
!->Remarks.-
!
!  1) LAPACK routines used:
!
!      -DGGEV:
!          Computes for a pair of N-by-N complex nonsymmetric matrices
!          (A,B), the generalized eigenvalues, and optionally, the left and/or
!          right generalized eigenvectors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine lapack_nonsymeigv (N,LenA,A,LenB,B,NumLambda,Lambda,Vectors,Info)
!!  end subroutine lapack_nonsymeigv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module interface_lapack
 
