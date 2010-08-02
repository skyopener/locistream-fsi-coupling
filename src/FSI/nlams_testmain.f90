module maindata
implicit none
!// Related to Mesh
       integer,save:: mainNumNodes
       real(8),save,dimension(:),allocatable::mainXgl,mainYgl,mainZgl

!// Related to Element
       integer,save:: mainNumElems
       integer,save,dimension(:,:),allocatable:: mainConnect

!// Related to BC
       integer,save:: mainNumBC
       integer,save,dimension(:),allocatable:: mainBCdof
       real(8),save,dimension(:),allocatable:: mainBCval

       integer,save:: mainAnsSize
       integer,save:: maintipNode 

!// Related to Material Property
!       real(8),save :: mainE1, mainNU12, mainStruDens,mainStruThickness
       real(8), parameter :: mainE1 = 70.d+9
       real(8), parameter :: mainNU12=0.3d0
       real(8), parameter :: mainStruDens= 2700.d0
       real(8), parameter :: mainStruThickness= 2.d-4

!// Related to Integration Scheme
       integer,save :: mainIntScheme
       real(8),parameter :: mainBeta = 0.25d0 , mainGanma= 0.5d0
       real(8),parameter :: mainGenAlph= 0.4d0

!// Related to displacement, velocity, accelaration
       real(8),save, dimension(:,:),allocatable :: maindisp
       real(8),save, dimension(:,:),allocatable :: mainvelo
       real(8),save, dimension(:,:),allocatable :: mainacce
       real(8),save, dimension(:,:),allocatable :: maindispRemesh

!// Related to aerodynamic force
       real(8),save, dimension(:,:),allocatable :: mainaeroforce
       real(8),save, dimension(:,:),allocatable :: mainaeroforce_pre

!// Related to inertial force
!       real(8),save, allocatable(:,1) :: mainforce_mbd

!// Related to internal force
!       real(8),save, allocatable(:,1) :: mainintforce_g

!// Related to Nodal_CSYS
       real(8),save, dimension(:,:,:,:),allocatable :: mainNodal_CSYS

!// Related to Time Advacing
       real(8),parameter :: maindeltaT= 5.d-6
       integer,save :: maintimeStep

!// Related to prescribed motion
       integer, save :: mainRotType, mainExtype, mainPlugType
       real(8), save :: mainFreq, mainPitchAmp, mainFlapAmp, mainLagAmp 
       real(8), save :: mainXfiAmp, mainYfiAmp, mainZfiAmp

      contains
       subroutine deallocatemaindata

       deallocate(mainXgl,mainYgl,mainZgl)
       deallocate(mainConnect)
       deallocate(mainBCdof,mainBCval)
       deallocate(maindisp,mainvelo,mainacce)
       deallocate(mainaeroforce,mainaeroforce_pre)
       deallocate(mainNodal_CSYS)
       deallocate(maindispRemesh)        
       end subroutine deallocatemaindata

end module maindata


