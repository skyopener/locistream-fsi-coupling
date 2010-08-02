      subroutine maininitilize
      use maindata

      implicit none
      mainAnsSize = mainNumNodes * 6
      maintipNode = 1

!      mainE1 = 70000000000.
!      mainE2 = mainE1
!      mainNU12 = 0.3
!      mainStruDens = 2700.
!      mainStruThickness = 0.0002
!      print *, 'thickness: ',mainStruThickness

      mainIntScheme = 2
!      mainBeta = 0.25d0
!      mainGanma= 0.50d0
!      mainGenAlph= 0.4d0

      allocate(maindisp(mainAnsSize,1))
      allocate(mainvelo(mainAnsSize,1))
      allocate(mainacce(mainAnsSize,1))

      maindisp = 0.0
      mainvelo = 0.0
      mainacce = 0.0

      allocate(mainaeroforce(mainAnsSize,1))
      allocate(mainaeroforce_pre(mainAnsSize,1))

      mainaeroforce =0.0
      mainaeroforce_pre =0.0

!      mainforce_mbd = 0.0

!      mainintforce_g = 0.0
      
      allocate(mainNodal_CSYS(3,mainNumElems,3,3));

!      mainNodal_CSYS = 
      mainNodal_CSYS(:,:,1,1) = 1.0
      mainNodal_CSYS(:,:,1,2) = 0.0
      mainNodal_CSYS(:,:,1,3) = 0.0

      mainNodal_CSYS(:,:,2,1) = 0.0
      mainNodal_CSYS(:,:,2,2) = 1.0
      mainNodal_CSYS(:,:,2,3) = 0.0

      mainNodal_CSYS(:,:,3,1) = 0.0
      mainNodal_CSYS(:,:,3,2) = 0.0
      mainNodal_CSYS(:,:,3,3) = 1.0

!      maindeltaT = 5.d-6 
      
!     mainExcite = 1
      mainExtype = 1 ! 0: plung, 1: flap

      mainFreq     = 30.0
! Flapping motion:
      mainRotType = 2 ! 1: sine, 2: 1-cos, 3: cos, 4: 1-exp*(-value*t^2))*sin
      mainPitchAmp =  0.0
      mainFlapAmp  = 17.0
      mainLagAmp   =  0.0
! Pluging motion:
      mainPlugType = 2 ! 1: 1-cos, 2: sin, 3: cos
      mainXfiAmp   =  0.0
      mainYfiAmp   =  0.0
      mainZfiAmp   =  0.0

      allocate(maindispRemesh(mainNumNodes,3)); maindispRemesh = 0.0
!      allocate(mainNLAMSdisp(mainNumNodes,3)); mainNLAMSdisp = 0.0;

      end subroutine maininitilize
