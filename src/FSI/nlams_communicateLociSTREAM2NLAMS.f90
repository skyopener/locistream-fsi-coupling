!        subroutine ComunicateLociSTREAM2NLAMS(A,B,C,D,E,F,G,H,I,J,K,L1,L2,M,N,O,P,R,S,T,U,V1,V2,Z,AA,BB,CC,DD,EE,FF,fwk1,fwk2,fwk3,pwk1,pwk2,pwk3,outdispRemesh)
!        subroutine ComunicateLociSTREAM2NLAMS(A,B,C,D,E,F,G,H,I,J,K,L1,L2,M,N,O,P,R,S,T,U,V1,V2,AA,BB,CC,DD,EE,FF,fwk1,fwk2,fwk3,pwk1,pwk2,pwk3,outdispRemesh)
        subroutine pass(rank,itfsi,A,B,C,D,E,F,G,H,I,J,K,L1,L2,M,N,O,P,R,S,T,U,V1,V2,Z,AA,BB1,tipnode,CC,EE,FF,fwk1,fwk2,fwk3,pwk1,pwk2,pwk3,outdispRemesh,outdisp,outvelo,outacce,outnodal)

! subroutine ComunicateLociSTREAM2NLAMS(A,B,C,D,E,F,J,K,L1,L2,M,N,O,P,R,S,T,U,V1,V2,W,Y,Z,AA,BB,CC,DD,EE,output)
!        subroutine ComunicateLociSTREAM2NLAMS(A, B)
        use LociSTREAMdata
        use nlams_excute
        use nlams_outdisp
        use nlams_init
!        use maindata

         implicit none
         integer:: i1
         integer:: csdfx,csdfy,csdfz,csdmx,csdmy,csdmz,tipnode
         integer:: rank, itfsi ! Loci-STREAM proc rank & itfsi
         integer:: A !// LociSTREAMNumNodes
         real(8) :: B(A),C(A),D(A) !// LocSTREAMXgl,LociSTREAMYgl,LociSTREAMZgl

         integer:: E !// LociSTREAMNumElems
         integer, dimension(E,3) :: F !// LociSTREAMConnect

         integer:: G !// LociSTREAMNumBC
         integer, dimension(G) :: H !// LociSTREAMBCdof
         real(8), dimension(G) :: I !// LociSTREAMBCval

         real(8) :: J, K, L1, L2 ! LociSTREAME1, LociSTREAMNU12, LociSTREAMStruDens, LociSTREAMStruThickness

         integer :: M ! LociSTREAMIntScheme
         real(8) :: N,O,P ! LociSTREAMBeta, LociSTREAMGanma, LociSTREAMGenAlph

         integer :: R ! LociSTREAMAnsSize
         real(8), dimension(R,1):: S, T, U !// LociSTREAMdisp, LociSTEAMvelo, LociSTEAMacce,
       
         real(8), dimension(R,1):: V1, V2  !// LociSTREAMaeroforce, LociSTREAMaeroforce_pre

         real(8), dimension(3,E,3,3) :: Z !// LociSTREAMNodal_CSYS

         real(8), dimension(R,1):: outdisp, outvelo, outacce
         real(8), dimension(A,3):: outdispRemesh
         real(8), dimension(3,E,3,3) :: outnodal

         real(8) :: AA !// DeltaT
         integer(8) :: BB1,BB2, CC, EE !// TimeStep, RotType, Extype, PlugType
         real(8):: FF ! Freq
         real(8):: fwk1,fwk2,fwk3,pwk1,pwk2,pwk3

!//
         character(14) :: SOLNAME

         LociSTREAMRank = rank ; LociSTREAMitfsi = 1
         if(LociSTREAMRank.eq.0) print *, 'rank & itfsi: ', LociSTREAMRank, LociSTREAMitfsi
         
         LociSTREAMPitchAmp = fwk1
         LociSTREAMFlapAmp  = fwk2
         LociSTREAMLagAmp   = fwk3

         LociSTREAMXfiAmp   = pwk1
         LociSTREAMYfiAmp   = pwk2
         LociSTREAMZfiAmp   = pwk3

         if (rank==0) print  *, "LociSTREAMPitchAmp, LociSTREAMFlapAmp, LociSTREAMLagAmp: ", LociSTREAMPitchAmp, LociSTREAMFlapAmp, LociSTREAMLagAmp
         if (rank==0) print  *, "LociSTREAMXfiAmp, LociSTREAMYfiAmp, LociSTREAMZfiAmp: ", LociSTREAMXfiAmp, LociSTREAMYfiAmp, LociSTREAMZfiAmp
!         read(*,*)

         LociSTREAMNumNodes = A
         if (rank==0) print  *, "LociSTREAMNumNodes: ", LociSTREAMNumNodes
         allocate(LociSTREAMXgl(LociSTREAMNumNodes)); LociSTREAMXgl=0.0
         allocate(LociSTREAMYgl(LociSTREAMNumNodes)); LociSTREAMYgl=0.0
         allocate(LociSTREAMZgl(LociSTREAMNumNodes)); LociSTREAMZgl=0.0

!         read(*,*)

!// Original
!         LociSTREAMXgl = B ; LociSTREAMYgl = C ; LociSTREAMZgl = D

!// pass values are defined w.r.t. CFD coordinate:
          LociSTREAMXgl = D ; LociSTREAMYgl = B ; LociSTREAMZgl = C
! 	if (rank==0) then
!          do i1=1, LociSTREAMNumNodes
!            print  *, LociSTREAMXgl(i1), LociSTREAMYgl(i1), LociSTREAMZgl(i1)
!          enddo
!          read(*,*)
!        endif

         LociSTREAMNumElems = E
         if (rank==0) print  *, "LociSTREAMNumElems: ", LociSTREAMNumElems
         allocate(LociSTREAMConnect(LociSTREAMNumElems,3))
         LociSTREAMConnect = F

!          if (rank==0) then
!          do i1=1, LociSTREAMNumElems
!            if (rank==0) print  *, 'LociConnect :',LociSTREAMConnect(i1,1:3)
!          enddo
!          read(*,*)
!         endif
 
         LociSTREAMNumBC = G
         if (rank==0) print  *, "LociSTREAMNumBC: ", LociSTREAMNumBC

         allocate(LociSTREAMBCdof(LociSTREAMNumBC))
         allocate(LociSTREAMBCval(LociSTREAMNumBC))
         
         LociSTREAMBCdof = H ; LociSTREAMBCval = I

!        if (rank==0) then
!          do i1=1, LociSTREAMNumBC
!            if (rank==0) print  *, 'LociSTREAMBCdof,LociSTREAMBCval :', LociSTREAMBCdof(i1),LociSTREAMBCval(i1)
!          enddo
!          read(*,*)  
!        endif

         LociSTREAME1 = J ; LociSTREAMNU12 = K ; LociSTREAMStruDens = L1 ; LociSTREAMStruThickness = L2

         if (rank==0) print  *, 'E1, NU12, Rho, Thickness :', LociSTREAME1,LociSTREAMNU12, LociSTREAMStruDens, LociSTREAMStruThickness
!         read(*,*)

         LociSTREAMIntScheme = M ; LociSTREAMBeta = N ; LociSTREAMGanma = O ; LociSTREAMGenAlph = P

         if (rank==0) print  *, 'IntScheme, Beta, Ganma, GenAlph :', LociSTREAMIntScheme,LociSTREAMBeta, LociSTREAMGanma, LociSTREAMGenAlph
!         read(*,*)

         LociSTREAMAnsSize = R
         if (rank==0) print  *, 'LociSTREAMAnsSize :', LociSTREAMAnsSize

         allocate(LociSTREAMdisp(LociSTREAMAnsSize,1)); LociSTREAMdisp= 0.0
         allocate(LociSTREAMvelo(LociSTREAMAnsSize,1)); LociSTREAMvelo= 0.0
         allocate(LociSTREAMacce(LociSTREAMAnsSize,1)); LociSTREAMacce= 0.0
         LociSTREAMdisp = S ; LociSTREAMvelo = T ; LociSTREAMacce = U

         allocate(LociSTREAMaeroforce(LociSTREAMAnsSize,1)); LociSTREAMaeroforce=0.0
         allocate(LociSTREAMaeroforce_pre(LociSTREAMAnsSize,1)); LociSTREAMaeroforce_pre= 0.0

         

!// V1 and V2 are defined w.r.t cfd coordinate:
         do i1=1,LociSTREAMNumNodes
         csdfx = (i1-1)*6+1
         csdfy = (i1-1)*6+2
         csdfz = (i1-1)*6+3
         csdmx = (i1-1)*6+4
         csdmy = (i1-1)*6+5
         csdmz = (i1-1)*6+6
      
         LociSTREAMaeroforce(csdfx,1) =  V1(csdfz,1)
         LociSTREAMaeroforce(csdfy,1) =  V1(csdfx,1)
         LociSTREAMaeroforce(csdfz,1) =  V1(csdfy,1)
         LociSTREAMaeroforce(csdmx,1) =  V1(csdmz,1)
         LociSTREAMaeroforce(csdmy,1) =  V1(csdmx,1)
         LociSTREAMaeroforce(csdmz,1) =  V1(csdmy,1)

         LociSTREAMaeroforce_pre(csdfx,1) =  V2(csdfz,1)
         LociSTREAMaeroforce_pre(csdfy,1) =  V2(csdfx,1)
         LociSTREAMaeroforce_pre(csdfz,1) =  V2(csdfy,1)
         LociSTREAMaeroforce_pre(csdmx,1) =  V2(csdmz,1)
         LociSTREAMaeroforce_pre(csdmy,1) =  V2(csdmx,1)
         LociSTREAMaeroforce_pre(csdmz,1) =  V2(csdmy,1)
         enddo

         

! TEMP: Check
!         LociSTREAMaeroforce = 0.0;
!         LociSTREAMaeroforce_pre = 0.0;

         if (rank==0) print  *, 'Writing out the forces:'
!         print *, "rank: ", rank
         
!// output aerodynamic force file
!// all variables are defined wrt csd coordinate:
  
         if (rank.eq.0) then
          open(10,file='locistreamaeroforce_wrt_csdcoordinate.dat')
          write(10,'(A)')'title = "Sample finite-element data"'
          write(10,*)'VARIABLES = "X", "Y", "Z", "FX", "FY", "FZ", "MX", "MY", "MZ","FXp", "FYp", "FZp", "MXp", "MYp", "MZp"'
          write(10,*)'ZONE F=FEPOINT, ET=TRIANGLE, N=',LociSTREAMNumNodes,',E=',LociSTREAMNumElems
          do i1=1,LociSTREAMNumNodes
           csdfx = (i1-1)*6+1
           csdfy = (i1-1)*6+2
           csdfz = (i1-1)*6+3
           csdmx = (i1-1)*6+4
           csdmy = (i1-1)*6+5
           csdmz = (i1-1)*6+6

          write (10,'(9(f18.9,1x))')LociSTREAMXgl(i1),LociSTREAMYgl(i1),LociSTREAMZgl(i1), LociSTREAMaeroforce(csdfx,1), LociSTREAMaeroforce(csdfy,1),  LociSTREAMaeroforce(csdfz,1) &
&  ,LociSTREAMaeroforce(csdmx,1), LociSTREAMaeroforce(csdmy,1), LociSTREAMaeroforce(csdmz,1)
         enddo
       
         do i1=1,LociSTREAMNumElems
          write (10,'(i6,1x,i6,1x,i6)')LociSTREAMConnect(i1,1),LociSTREAMConnect(i1,2),LociSTREAMConnect(i1,3)
         enddo
        close(10)
       endif

       if (rank.eq.0) then
         open(10,file='locistreamdispvelacce.dat')
         do i1=1,LociSTREAMAnsSize
         write(10,*) LociSTREAMdisp(i1,1),LociSTREAMvelo(i1,1),LociSTREAMacce(i1,1)
         enddo
         close(10)
       endif
         
!         if (rank==0) print  *, 'After writing outfiles:'
                 
         allocate(LociSTREAMNodal_CSYS(3,LociSTREAMNumElems,3,3)) ; LociSTREAMNodal_CSYS=0.0
         LociSTREAMNodal_CSYS = Z

!         if (rank==0) print  *, LociSTREAMNodal_CSYS(1,1,1,1), LociSTREAMNodal_CSYS(1,1,1,2),LociSTREAMNodal_CSYS(1,1,1,3)
!         if (rank==0) print  *, LociSTREAMNodal_CSYS(1,1,2,1), LociSTREAMNodal_CSYS(1,1,2,2),LociSTREAMNodal_CSYS(1,1,2,3)
!         if (rank==0) print  *, LociSTREAMNodal_CSYS(1,1,3,1), LociSTREAMNodal_CSYS(1,1,3,2),LociSTREAMNodal_CSYS(1,1,3,3)
!         read(*,*)

         LociSTREAMdeltaT = AA; LociSTREAMtimeStep = BB1 !; LociSTREAMtimeStepfsi = BB2
         LociSTREAMRotType = CC !; LociSTREAMExtype = DD
         LociSTREAMPlugType = EE ; LociSTREAMFreq = FF

!         if (rank==0) print  *, 'DeltaT, timeStep,Rotype, Extype,PlugType,Freq :', &
!&LociSTREAMdeltaT, LociSTREAMtimeStep, LociSTREAMRotType, LociSTREAMExtype, LociSTREAMPlugType, LociSTREAMFreq
         if (rank==0) print  *, 'DeltaT, timeStep ,Rotype, PlugType, Freq :', &
&LociSTREAMdeltaT, LociSTREAMtimeStep, LociSTREAMRotType, LociSTREAMPlugType, LociSTREAMFreq

         LociSTREAMtipNode = tipnode !18 for square2d 
         if (rank==0) print  *, 'tipnode;',tipnode

!         if(LociSTREAMtimeStep.eq.1)then 
!         if(LociSTREAMtimeStepfsi.eq.1)then
         allocate(NLAMSdisp(LociSTREAMAnsSize,1)); NLAMSdisp = 0.0d0
         allocate(NLAMSvelo(LociSTREAMAnsSize,1)); NLAMSvelo = 0.0d0
         allocate(NLAMSacce(LociSTREAMAnsSize,1)); NLAMSacce = 0.0d0
         allocate(NLAMSNodal_CSYS(3,LociSTREAMNumElems,3,3)); NLAMSNodal_CSYS= 0.0d0
         allocate(NLAMSdispRemesh(LociSTREAMNumNodes,3)); NLAMSdispRemesh = 0.0d0
!         endif

         call ExcuteNLAMS
         if (rank==0) print  *, 'after ExcuteNLAMS'
!         call NLAMStoLociSTREAM 

         outdisp=0.0d0; outvelo=0.0d0; outacce=0.0d0; outnodal=0.0d0; outdispRemesh=0.0d0
         outdisp = NLAMSdisp
         outvelo = NLAMSvelo
         outacce = NLAMSacce
         outnodal = NLAMSNodal_CSYS
         
! // outdispRemesh is defined w.r.t cfd coordinate:
         outdispRemesh = NLAMSdispRemesh
   
        if (rank.eq.0) then
         WRITE(SOLNAME,"('disp',I6.6,'.dat')"),LociSTREAMtimeStep
         open(10,file=SOLNAME)
        do i1=1,LociSTREAMNumNodes
          write(10,'(i6,1x,3(e14.7,1x))') i1,outdispRemesh(i1,1),outdispRemesh(i1,2),outdispRemesh(i1,3)
        enddo
         close(10)
       endif

!         if (rank==0) print  *, "in communicate"
!         if (rank==0) print  *, NLAMSNodal_CSYS(1,1,1,1),NLAMSNodal_CSYS(1,1,1,2),NLAMSNodal_CSYS(1,1,1,3)
!         if (rank==0) print  *, NLAMSNodal_CSYS(1,1,2,1),NLAMSNodal_CSYS(1,1,2,2),NLAMSNodal_CSYS(1,1,2,3)
!         if (rank==0) print  *, NLAMSNodal_CSYS(1,1,3,1),NLAMSNodal_CSYS(1,1,3,2),NLAMSNodal_CSYS(1,1,3,3)
!         if (rank==0) print  *, outnodal(1,1,1,1),outnodal(1,1,1,2),outnodal(1,1,1,3)
!         if (rank==0) print  *, outnodal(1,1,2,1),outnodal(1,1,2,2),outnodal(1,1,2,3)
!         if (rank==0) print  *, outnodal(1,1,3,1),outnodal(1,1,3,2),outnodal(1,1,3,3)


         call LociSTREAMdata_deallocation
         if (rank==0) print  *, 'after LociSTREAMdata_deallocation'


! end subroutine ComunicateLociSTREAM2NLAMS
end subroutine pass
