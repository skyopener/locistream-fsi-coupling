module LociSTREAMdata
implicit none
!// Related to Loci-STREAM rank
       integer,save:: LociSTREAMRank, LociSTREAMitfsi

!// Related to Mesh
       integer,save:: LociSTREAMNumNodes
       real(8),save,dimension(:),allocatable::LociSTREAMXgl,LociSTREAMYgl,LociSTREAMZgl

!// Related to Element
       integer,save:: LociSTREAMNumElems
       integer,save,dimension(:,:),allocatable:: LociSTREAMConnect

!// Related to BC
       integer,save:: LociSTREAMNumBC
       integer,save,dimension(:),allocatable:: LociSTREAMBCdof
       real(8),save,dimension(:),allocatable:: LociSTREAMBCval

       integer,save:: LociSTREAMAnsSize
       integer,save:: LociSTREAMtipNode 

!// Related to Material Property
       real(8),save :: LociSTREAME1, LociSTREAMNU12, LociSTREAMStruDens,LociSTREAMStruThickness

!// Related to Integration Scheme
       integer,save :: LociSTREAMIntScheme
       real(8),save :: LociSTREAMBeta, LociSTREAMGanma
       real(8),save :: LociSTREAMGenAlph

!// Related to displacement, velocity, accelaration
       real(8),save, dimension(:,:),allocatable :: LociSTREAMdisp
       real(8),save, dimension(:,:),allocatable :: LociSTREAMvelo
       real(8),save, dimension(:,:),allocatable :: LociSTREAMacce
!       real(8),save, allocatable :: LociSTREAMdispRemesh(:,:)

!// Related to aerodynamic force
       real(8),save, dimension(:,:),allocatable :: LociSTREAMaeroforce
       real(8),save, dimension(:,:),allocatable :: LociSTREAMaeroforce_pre

!// Related to inertial force
!       real(8),save, allocatable(:,1) :: LociSTREAMforce_mbd

!// Related to internal force
!       real(8),save, allocatable(:,1) :: LociSTREAMintforce_g

!// Related to Nodal_CSYS
       real(8),save, dimension(:,:,:,:),allocatable :: LociSTREAMNodal_CSYS

!// Related to Time Advacing
       real(8),save :: LociSTREAMdeltaT
       integer,save :: LociSTREAMtimeStep

!// Related to prescribed motion
       integer, save :: LociSTREAMRotType, LociSTREAMExtype, LociSTREAMPlugType
       real(8), save :: LociSTREAMFreq, LociSTREAMPitchAmp, LociSTREAMFlapAmp, LociSTREAMLagAmp 
       real(8), save :: LociSTREAMXfiAmp, LociSTREAMYfiAmp, LociSTREAMZfiAmp

!// Related to Output
       real(8),save, dimension(:,:),allocatable :: NLAMSdisp
       real(8),save, dimension(:,:),allocatable :: NLAMSvelo
       real(8),save, dimension(:,:),allocatable :: NLAMSacce
       real(8),save, dimension(:,:,:,:),allocatable :: NLAMSNodal_CSYS
       real(8),save, dimension(:,:), allocatable :: NLAMSdispRemesh

      contains
        subroutine LociSTREAMdata_deallocation
         deallocate(NLAMSdisp,NLAMSvelo,NLAMSacce,NLAMSNodal_CSYS,NLAMSdispRemesh)
         deallocate(LociSTREAMXgl,LociSTREAMYgl,LociSTREAMZgl)
         deallocate(LociSTREAMConnect)
         deallocate(LociSTREAMBCdof, LociSTREAMBCval)
         deallocate(LociSTREAMdisp, LociSTREAMvelo, LociSTREAMacce)
         deallocate(LociSTREAMaeroforce, LociSTREAMaeroforce_pre)
         deallocate(LociSTREAMNodal_CSYS)
        end subroutine LociSTREAMdata_deallocation
     end module LociSTREAMdata
