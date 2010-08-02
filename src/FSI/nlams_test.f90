      program main
      use maindata

      implicit none
!      use maindata
      integer :: timestep
      integer :: i

      call readgrid
      print *, 'after main_readgrid'
      call readBC
      print *, 'after main_readBC'
      call maininitilize
      print *, 'after main_maininitilize'

      do timestep=1,7000
! Comunication between Loci-STREAM to NLAMS
! Then run NLAMS
      print *, 'timestep: ',timestep

!      open(10,file='grid_ori.dat')
!      do i=1,mainNumNodes
!      write(10,*), mainXgl(i),mainYgl(i),mainZgl(i)
!      enddo
!      close(10)

      call ComunicateLociSTREAM2NLAMS(mainNumNodes,mainXgl,mainYgl,mainZgl,mainNumElems,mainConnect,mainNumBC,mainBCdof,mainBCval,mainE1,mainNU12,mainStruDens,mainStruThickness,mainIntScheme,mainBeta,mainGanma,mainGenAlph,mainAnsSize,maindisp,mainvelo,mainacce,mainaeroforce,mainaeroforce_pre,mainNodal_CSYS,maindeltaT,timestep,mainRotType,mainExtype,mainPlugType,mainFreq,mainPitchAmp,mainFlapAmp,mainLagAmp,mainXfiAmp,mainYfiAmp,mainZfiAmp,maindispRemesh)

      if(timestep.eq.1)then 
      open(10,file='acce_timestep1.dat')
      do i=1,mainAnsSize
      write(10,*) mainacce(i,1)
      enddo
      close(10)
      elseif(timestep.eq.2)then
      open(10,file='disp_timestep2.dat')
      do i=1,mainAnsSize
      write(10,*) maindisp(i,1)
      enddo
      close(10)
      elseif(timestep.eq.3)then
      open(10,file='disp_timestep3.dat')
      do i=1,mainAnsSize
      write(10,*) maindisp(i,1)
      enddo
      close(10)
      endif
      enddo

      call deallocatemaindata
      print *, 'after maindata_deallocation'

      stop
      end program main


