      module nlams_outdisp

      contains
       subroutine NLAMStoLociSTREAM
         use mesh_connect_vars, only: NumNodes, UndefGrid
         use nlams_init, only: FinalAns,AnsSize,D_curr,Ddot_curr,Dddot_curr,Nodal_CSYS
         use nlams_input, only: PlugType, RotType
         use LociSTREAMdata, only: NLAMSdisp,NLAMSvelo,NLAMSacce,NLAMSNodal_CSYS,NLAMSdispRemesh,LociSTREAMRank
         use nlams_undeformed,   only: MemOri, T
         
         implicit none

!         real(8), dimension(AnsSize,1),intent(in):: a
!         real(8), dimension(NumNodes,3),intent(out):: b
!         real(8):: a(AnsSize,1)

         integer:: i,jstart,jend, j1, j2, j3
         integer:: kstart, kend
         real(8):: TotDisp(NumNodes,3), CurrCoord(NumNodes,3), CurrCoordInert(NumNodes,3), Aux(3,1)
         real(8):: Tdisp(AnsSize,1) 

!         allocate(NLAMSdisp(NumNodes,3))
! CFD: x is chordwise, y is vertical, and z is spanwise
! CSD: x is spanwise, y is chordwise, and z is normal to plane
! so transform here:
          TotDisp = 0.0d0
          kstart = 1; kend = 3

!          do i=1,NumNodes
!           TotDisp(i,1:3) = Tdisp(kstart:kend,1)
!           kstart = 6*(i-1)+ 1
!           kend   = kstart + 2
!          enddo

!          if(LociSTREAMRank==0) then
!          do i=1,3
!          print *, T(i,1:3)
!          enddo
!          endif

!          if(PlugType.ne.0.and.RotType.eq.0)then
!          T      = 0.d0
!          T(1,1) = 1.d0
!          T(2,2) = 1.d0
!          T(3,3) = 1.d0
!          endif

            CurrCoord = 0.d0 ; CurrCoordInert = 0.d0

!          do i=1,NumNodes
!            j1 = (i-1)*6 + 1
!            j2 = (i-1)*6 + 2
!            j3 = (i-1)*6 + 3
!	    CurrCoord(i,1)     = TotDisp(i,j1)
!	    CurrCoord(i,2)     = TotDisp(i,j2)
!	    CurrCoord(i,3)     = TotDisp(i,j3)
!		   Aux(:,1)     = CurrCoord(i,1:3)
!                   Aux(:,1)     = matmul(T,Aux(:,1))
!           CurrCoordInert(i,1:3)= Aux(1:3,1)
!            NLAMSdispRemesh(i,1) = CurrCoordInert(i,2) + MemOri(2,1)
!            NLAMSdispRemesh(i,2) = CurrCoordInert(i,3) + MemOri(3,1)
!            NLAMSdispRemesh(i,3) = CurrCoordInert(i,1) + MemOri(1,1)
!          enddo



!           print *, 'ABC'

!           do i=1,AnsSize
!           if(FinalAns(i).le.1.d-18)then
!            FinalAns(i) = 0.d0
!           endif
!           enddo


         do i=1,NumNodes
          j1 = (i-1)*6 + 1
          j2 = (i-1)*6 + 2
          j3 = (i-1)*6 + 3
           NLAMSdispRemesh(i,1) = FinalAns(j2) + MemOri(2,1)
           NLAMSdispRemesh(i,2) = FinalAns(j3) + MemOri(3,1)
           NLAMSdispRemesh(i,3) = FinalAns(j1) + MemOri(1,1)

! Testing rigid body motion only for fsi framework
!          NLAMSdispRemesh(i,1) = MemOri(2,1)
!          NLAMSdispRemesh(i,2) = MemOri(3,1)
!          NLAMSdispRemesh(i,3) = MemOri(1,1)
         enddo
!          if(LociSTREAMRank==0) print *, 'B2', LociSTREAMRank
!         print *,  'before nlams2locistreamdispvelacce written'

! For CSD simulation at next time step         
         NLAMSdisp = D_curr
!          if(LociSTREAMRank==0) print *, 'B3', LociSTREAMRank
         NLAMSvelo = Ddot_curr
!          if(LociSTREAMRank==0) print *, 'B4', LociSTREAMRank
         NLAMSacce = Dddot_curr
!          if(LociSTREAMRank==0) print *, 'B5', LociSTREAMRank
         NLAMSNodal_CSYS = Nodal_CSYS
!          if(LociSTREAMRank==0) print *, 'B6', LociSTREAMRank
!         print *,  'right before nlams2locistreamdispvelacce written'

         if(LociSTREAMRank==0)then
          open(10,file='nlams2locistreamdisp.dat')
          open(20,file='nlams2locistreamvelo.dat')
          open(30,file='nlams2locistreamacce.dat')
!          write(10,*) "Displacement, Velocity, Acceralation"

          do i=1,NumNodes
           jstart = 6*(i-1)+1
           jend   = jstart + 5
           write(10,'(i4,6(1x,f14.7))') i, NLAMSdisp(jstart:jend,1)
           write(20,'(i4,6(1x,f14.7))') i, NLAMSvelo(jstart:jend,1)
           write(30,'(i4,6(1x,f14.7))') i, NLAMSacce(jstart:jend,1)
          enddo
          close(10); close(20); close(30)
         endif

!         if(LociSTREAMRank==0)print *, 'B8', LociSTREAMRank
!       print *,  'nlams2locistreamdispvelacce written'
         
       end subroutine NLAMStoLociSTREAM
      end module nlams_outdisp
