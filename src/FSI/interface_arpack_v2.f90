!->Copyright by The University of Michigan, Aerospace Department. 2006
!
!->Module INTERFACE_ARPACK. Rafa Palacios. 04Oct2006
!
!->Description.-
!
!  This module includes interface for ARPACK routines.
!
!->Subroutines:
!
!
!->Remarks.-
!
!  1) External routines are part of the LAPACK 3 library, which needs
!     to be compiled with this module.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module interface_arpack
 use mod_shared
 implicit none

!   (There are no public variables in this module).

 contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ARPACK_NONSYMEIGV
!
!->Description.-
!
!   Computes the complex eigenvalues of the generalized problem AX= lambda*B*X.
!
!->Remarks.-
!
!  1) ARPACK routines used:
!
!      -dnaupd: Implicitly Restarted Arnoldi iteration.
!      -dneupd: Postprocess eigenvalues and eigenvectors.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine arpack_nonsymeigv (N,LenA,A,LenB,B,NumLambda,Lambda,Vectors,Info)
  use lib_sparse
  use interface_lapack

  integer,     intent(in) :: N           ! Size of the matrix.
  integer,     intent(in) :: LenA        ! Size of the storage array for A matrix.
  type(sparse),intent(in) :: A (:)       ! Matrix A.
  integer,     intent(in) :: LenB        ! Size of the storage array for B matrix.
  type(sparse),intent(in) :: B (:)       ! Matrix B, definite positive.
  integer,     intent(in) :: NumLambda   ! Number of output eigenvalues.
  complex(8),  intent(out):: Lambda(:)   ! Complex Eigenvalues.
  complex(8),  intent(out):: Vectors(:,:)! Eigenvectors.
  integer,     intent(out):: Info        ! Error codes.

! Local variables.

  real(8),allocatable:: ArnV(:,:)       ! Arnoldi basis vectors.
  real(8),allocatable:: Bband(:,:)      ! B matrix in bandwidth form.
  real(8),allocatable:: D(:,:)
  integer            :: DimArnV         ! Dimension of Arnoldi vector space.
  integer            :: DimWorkL        ! Dimension of Working Arrary WorkL.
  integer            :: EigPoint(14)    ! Pointer array for eigenv calculation.
  integer            :: EigParams(11)   ! Input parameters for eigv calculation.
  real(8),allocatable:: EigResid(:)     ! Residual vector in eigv calculation.
  logical,allocatable:: EigSelect(:)    ! Not used.
  real(8)            :: EigTol          ! Tolerance in Eigenvalue calculation.
  real(8),allocatable:: EigVec(:,:)     ! Eigenvector matrix given by dneupd.
  logical            :: Flag            ! Logical flag :).
  integer            :: i,j,k           ! Counters.
  complex(8)         :: iu              ! Imaginary Unit.
  integer            :: ido             ! Reverse communication flag.
  integer            :: KL              ! Number of subdiagonals in the band of B.
  integer            :: KU              ! Number of superdiagonals in the band of B.
  integer            :: LenB1A          ! Used length in storage array for B1A matrix.
  integer,allocatable:: List(:)         ! Ordered list of eigenvalues.
  integer,allocatable:: LUPoint(:)      ! Pointer array for LU decomposition.
  real(8)            :: MaxIm           ! Maximum of imaginary values.
  real(8)            :: SigmaR,SigmaI
  real(8),allocatable:: WorkD(:)        ! Working Arrays.
  real(8),allocatable:: WorkE(:)     
  real(8),allocatable:: WorkL(:)

  iu= cdsqrt(dcmplx(-1.D0))


! Compute LU decomposition of B.
!==============================

  allocate(LUPoint(N))

! Determine upper (KU) and lower (KL) bandwidths of matrix B.

  call sparse_getbandwidth (LenB,B,KL,KU)

! Convert input matrix B to banded form and get LU decomposition.

  allocate(Bband(2*KL+KU+1,N)); Bband  = 0.d0
  call sparse_copy2band (LenB,B,KL,KU,BBand)
  call lapack_band_lufact (N,KL,KU,Bband,LUPoint,Info)
  if (Info.ne.0) return

! Compute Eigenvalues of largest imaginary part.
!==============================================

! Set up workspace.

  DimArnV = NumLambda+4
  DimWorkL= 3*DimArnV**2+6*DimArnV

  allocate(EigResid(N));        EigResid = 0.d0
  allocate(ArnV(N,DimArnV));    ArnV     = 0.d0
  allocate(WorkD(3*N));         WorkD    = 0.d0
  allocate(WorkL(DimWorkL));    WorkL    = 0.d0

! Iterative estimation of the eigenvalues (ARPACK).
! Options:
! * 'I': Standard eigenvalue problem A*x = lambda*x
! * 'LI': Compute eigenvalues of largest imaginary part.

  ido = 0
  info= 0
  EigTol = 0.d0
  EigParams(1)= 1
  EigParams(3)= 300   ! Max Number of iterations.
  EigParams(7)= 1

  Flag=.true.
  do while (Flag.or.(ido.eq.-1).or.(ido.eq.1))
    Flag=.false.
    call dnaupd (ido,'I',N,'LI',NumLambda,EigTol,EigResid, &
&                DimArnV,ArnV,N,EigParams,EigPoint,WorkD,  &
&                WorkL,DimWorkL,Info)

    call arpack_nonsymeigv_fun (N,LenA,A,KL,KU,Bband,LUPoint,&
&                               WorkD(EigPoint(1):EigPoint(1)+N-1), &
&                               WorkD(EigPoint(2):EigPoint(2)+N-1))

  end do

  if (Info.ne.0) return
  deallocate (LUPoint,Bband)

! Postprocess eigenvalues and eigenvectors.
! Options:
! * Flag=.true.: Eigenvectors are computed.
! * 'A': Compute Ritz vectors.

  allocate(D(NumLambda+1,2));      D        = 0.d0
  allocate(EigSelect(DimArnV));    EigSelect=.true.
  allocate(EigVec(N,NumLambda+1)); EigVec   = 0.d0
  allocate(WorkE(3*N));            WorkE    = 0.d0

  call dneupd (.true.,'A',EigSelect,D(:,1),D(:,2),EigVec,N,SigmaR,  &
&              SigmaI,WorkE,'I',N,'LI',NumLambda,EigTol,   & 
&              EigResid,DimArnV,ArnV,N,EigParams,EigPoint, &
&              WorkD,WorkL,DimWorkL,Info)

  if (Info.ne.0) return
  deallocate (EigSelect,EigResid,WorkD,WorkE,WorkL)

! Identify the larger M eigenvalues.

  allocate(List(NumLambda)); List  =0

  Vectors=(0.d0,0.d0)
  do i=1,NumLambda
    MaxIm=0.d0

    do j=1,NumLambda
      Flag=.true.
      do k=1,i-1
        if (List(k).eq.j) Flag=.false.
      end do

      if (Flag) then
        if (abs(D(j,2)).gt.MaxIm) then
          MaxIm=abs(D(j,2))
          List(i)=j
        end if
      end if
    end do

! Store eigenvalues.

    Lambda(i)=   dcmplx(D(List(i),1)) + iu*dcmplx(D(List(i),2))
  end do

! Store eigenvectors.

  i=1
  do while (i.le.NumLambda)
    Vectors(:,i)  =dcmplx(EigVec(:,List(i))) + iu*dcmplx(EigVec(:,List(i+1)))
    Vectors(:,i+1)=dcmplx(EigVec(:,List(i))) - iu*dcmplx(EigVec(:,List(i+1)))
    i=i+2
  end do

! End of routine.

  deallocate (D,List,EigVec)

  return
 end subroutine arpack_nonsymeigv



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine arpack_nonsymeigv_fun (N,LenA,A,KL,KU,Bband,iPoint,v,w)
  use lib_sparse

  integer,     intent(in):: N           ! Dimension of the problem.
  integer,     intent(in):: LenA        ! Size of the storage array for matrix A.
  type(sparse),intent(in):: A (:)       ! Matrix A.
  integer,     intent(in):: KL,KU       ! Bandwidth dimensions of matrix B.
  real(8),     intent(in):: Bband(:,:)  ! Matrix B (banded form).
  integer,     intent(in):: iPoint(:)   ! Pointer array.
  real(8),     intent(in):: v(:)        ! Input vector.
  real(8),     intent(out)::w(:)        ! Output vector.
  
  integer::i,j,k                        ! Counters.
  integer::iErr                         ! Error Code
  real(8),allocatable:: RHS(:)          ! RHS for the inversion process.
 
! Initialize.

  allocate(RHS(N))
  w=0.d0

! Loop on the columns of matrix A.

  do j=1,N
    RHS=0.d0
    do k=1,LenA
      if (A(k)%j.eq.j) RHS(A(k)%i)=A(k)%a
    end do

    call DGBTRS ('N',N,KL,KU,1,Bband,2*KL+KU+1,iPoint,RHS,N,iErr)
    if (iErr.ne.0) call error (iuAllOut,ErrEigSol,'Error in LU backsubstitution')

    do i=1,N
      w(i)=w(i) + RHS(i)*v(j)
    end do
  end do

  deallocate(RHS)
 end subroutine arpack_nonsymeigv_fun
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module interface_arpack
