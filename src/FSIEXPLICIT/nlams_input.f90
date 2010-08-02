!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008

!-> Module.- NLAMS_INPUT. Satish Chimakurthi
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!    Get input data for NLAMS.
!
!->Subroutines:
!
!     |-nlams_input_soltype:   Type of solution.
!
!
!->Remarks:
!
!   1) The variable DataBlock get the information of the fields 2 to 9 (one
!      or several lines) in the input cards. When used in the get_next_entry
!      subroutine must have a length multiple of 64 (the multiple is number
!      of expected lines).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nlams_input
  use LociSTREAMdata
!  use mod_shared
  
  implicit none

  !Isotropic/Orthotropic material properties
!  integer,save             :: FlagIso
  integer, parameter :: FlagIso = 0 ! Isotropic

  real(8),save             :: E1, E2, NU12, NU21, G12, RHO
  
  !Layup information
  integer, parameter :: NumLayers = 0

!  integer,save             :: NumLayers
  real(8),save,dimension(8):: PlyAngle
  real(8),save             :: Thickness
  
  ! beta parameters for the OPT element formulation
  real(8),save:: alfab,beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8
  real(8),save:: beta9
  
  !vParam GAUSS card
  integer,save    :: NumGaussPoints

  !vParam ELEMENT card
  integer,save    :: EleFlag
  
  !vParam DEBUG card
  integer,save    :: DebugPar
 
!
  integer,save      :: Iter,Llevel
  integer,parameter :: MaxLlevel = 1 !Maximum number of load levels
  integer,parameter :: MaxIter = 10000  !Maximum number of NR iterations
  integer,parameter :: OutputIter = 200 !Output tecplot/ascii files after this iteration
  integer,parameter :: NumTSteps = 150000  !Maximum number of time-steps

  !integer,save :: StressStiff  !In the input file now
  integer, parameter :: StressStiff = 1
  real(8),parameter :: Tol = 0.0005d0
  integer,save      :: ConvFlag = 0
  real(8),save      :: PrevNorm

  !Dynamics
  integer,save ::TStep,NSubIter 
  real(8),save ::num

  !Damping constant
  real(8),save   ::DampCons = 0.0d0

  !sine/cos excitation type
  integer,save   ::RotType  !If 1, a sine excitation will be applied else 1-cos  
  integer,save   ::ExType
  integer,save   ::PlugType

  !vParam TIPNODE
  integer,save   ::TipNode   !Output node at which displacement in the vertical direction will be written

  !vParam INTSCHEM
  integer,save   ::IntScheme  !1-> Newmark, 2->Gen.Alfa, 3->Houbolt

  !vParam TIMESTEP -> Time-step
  real(8),save :: DeltaT        !Timestep   !Now from the input file

  !For Newmark
  real(8),save :: GammaNewmark
  real(8),save :: BetaNewmark

!For single-step houbolt
!suggested: -0.5, 1.5
!real(8),save  :: gammahoub
!real(8),save  :: gammaphoub

  !For generalized-alpha
  real(8),save  :: SpRadius

!For wilson-theta
!real(8),save :: wilsontheta

  !For plunge motion
  real(8),save  :: XfiAmp,YfiAmp,ZfiAmp
  real(8),save  :: PlungeOmega !This value is in rad/sec

  !For flapping motion
  real(8),save    ::PitchAmp,FlapAmp,LagAmp
  real(8),save    ::FlapOmega

  !Rayleigh damping coefficients
  real(8),save::AlfaR     
  real(8),save::BetaR

  !Delta for convergence within the NR loop in the structural solver
  real(8),save:: CSDMinDelta  !Input file

contains
 
 subroutine NLAMS_input_LociSTREAM
  use LociSTREAMdata
!  use nlams_input

!  print *, "first in input"

  CSDMinDelta = 1e-6
  AlfaR       = 0.0
  BetaR       = 0.0

  FlapOmega = LociSTREAMFreq
  PitchAmp  = LociSTREAMPitchAmp ; FlapAmp  = LociSTREAMFlapAmp ; LagAmp  = LociSTREAMLagAmp

  PlungeOmega = LociSTREAMFreq
  XfiAmp    = LociSTREAMXfiAmp ; YfiAmp = LociSTREAMYfiAmp ; ZfiAmp = LociSTREAMZfiAMp

  TStep = LociSTREAMtimeStep
  DeltaT = LociSTREAMdeltaT
  IntScheme = LociSTREAMintScheme

  TipNode = LociSTREAMtipNode

  RotType = LociSTREAMRotType; ExType = 1 ; PlugType = LociSTREAMPlugType

  EleFlag = 1; NumGaussPoints = 3

! print *, "middle in input"

  alfab = 1.5d0; beta1 = 1.d0; beta2 = 2.d0; beta3 = 1.d0; beta4 = 0.0
  beta5 = 1.d0; beta6 = -1.d0; beta7 = -1.d0; beta8 = -1.d0; beta9 = -2.0
  
  SpRadius  = LociSTREAMGenAlph

  GammaNewmark = LociSTREAMGanma
  BetaNewmark  = LociSTREAMBeta
 
  PlyAngle = 0.0

  E1   = LociSTREAME1 ; E2 = E1 ;
!  E1 = 70d+9
  NU12 = LociSTREAMNU12 ; NU21  = LociSTREAMNU12
  G12  = 0.0
  RHO  = LociSTREAMStruDens
  Thickness = LociSTREAMStruThickness 

!  Thickness = 2d-4
! print *, "last in input"
 end subroutine NLAMS_input_LociSTREAM
end module nlams_input
