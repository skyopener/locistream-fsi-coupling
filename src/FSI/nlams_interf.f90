!-> Copyright by The University of Michigan, Aerospace Department. 2009
!
!-> Subroutine.- nlams_output
!
!-> Remarks.-
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_interf
! use mod_shared
 implicit none
!
!-> (There are no public variables in this module).
!
 contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Subroutine nlams_INTERF_TECPLOT
!
!-> Description:
!
!    Write a data on an ordered grid in TECPLOT ASCII format.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nlams_output (Tdisp)
use mesh_connect_vars, only: NumDofPerNode,NumNodes,Connect,UndefGrid,NumElems
use nlams_init,         only: AnsSize
use nlams_input,        only: Iter,OutputIter,TStep,TipNode
use nlams_undeformed,   only: T,MemOri
use shared_time,        only: TheTime
use LociSTREAMdata, only: LociSTREAMRank, LociSTREAMaeroforce,LociSTREAMitfsi

implicit none
!Local variables
integer :: kstart,kend,i, tipdof1, tipdof2, tipdof3
real(8),dimension(NumNodes,3)::TotDisp
real(8),dimension(NumNodes,3)::CurrCoord,CurrCoordInert
real(8),dimension(3,1)::Aux,Tdisp1

!intent(in) variables
real(8),dimension(AnsSize,1),intent(in)::Tdisp
!integer,intent(in)			  

!For tecplot output files
character(len=100)::outputFileNLAMS,casename
integer::unitnum

! Initialization
  tipdof1 = 0; tipdof2 = 0; tipdof3 = 0
  CurrCoord = 0.0d0
  CurrCoordInert =0.0d0

! TipNode is number of node givin in input file
! Tdisp is pure deformation vector with respect to global frame 
     tipdof1=((TipNode-1)*6)+1
     tipdof2=((TipNode-1)*6)+2
     tipdof3=((TipNode-1)*6)+3
!     Aux(1,1) = Tdisp(tipdof1,1)
!     Aux(2,1) = Tdisp(tipdof2,1)
!     Aux(3,1) = Tdisp(tipdof3,1)
!     Aux(1:3,1) = Aux(1:3,1) + UndefGrid(TipNode,1:3) + MemOri(1:3,1)
!
!     Tdisp1(:,1) = matmul(T,Aux(:,1))

!First col time, second col: displacement in the vertical dirn at the tip
!     write(204,*)TheTime,Tdisp1(3,1)

!Write energy
!     write(339,'(E14.6,1x,3E14.6)')TheTime,KinetE(1,1),PotenE,TotE

!Write ASCII output file of the current X,Y,Z coordinates w.r.t. the inertial frame
     TotDisp = 0.0d0
     kstart = 1
     kend = 3
     do i=1,NumNodes
      TotDisp(i,1:3) = Tdisp(kstart:kend,1)
      kstart = kstart + NumDofPerNode
      kend   = kend   + NumDofPerNode
     enddo
! if (LociSTREAMRank==0)  print *,   'A1', LociSTREAMitfsi
 
! Write origin of global coordinate system with respect to inertial frame
! If wing is pluging, MemOri(1:3,1) will be changed.
     if (LociSTREAMRank==0)  write(393,'(E14.6,1x,3E14.6)')TheTime,MemOri(1,1),MemOri(2,1),MemOri(3,1)
! if (LociSTREAMRank==0)  print *,   'A2', LociSTREAMitfsi
! Current coordinate system in global system
     do i=1,NumNodes
      CurrCoord(i,1:3)     = UndefGrid(i,1:3) + TotDisp(i,1:3) + MemOri(1:3,1)
              Aux(:,1)     = CurrCoord(i,1:3)
              Aux(:,1)     = matmul(T,Aux(:,1))
      CurrCoordInert(i,1:3)=Aux(1:3,1)
     enddo
! if (LociSTREAMRank==0)  print *,   'A3', LociSTREAMitfsi
!     if (LociSTREAMRank==0)  write(204,'(e14.7,1x,9e15.7)')TheTime,(CurrCoord(TipNode,i),i=1,3),(CurrCoordInert(TipNode,i),i=1,3),(TotDisp(TipNode,i),i=1,3)
if (LociSTREAMRank==0)  write(204,'(e14.7,1x,9e15.7)')TheTime,(CurrCoord(TipNode,i),i=1,3),(CurrCoordInert(TipNode,i),i=1,3), LociSTREAMaeroforce(tipdof1,1) &
&, LociSTREAMaeroforce(tipdof2,1),LociSTREAMaeroforce(tipdof3,1)

!if (LociSTREAMRank==0)  print *,   'A4', LociSTREAMitfsi

if (LociSTREAMRank==0) then
! print *, 'Tstep in output1 ', LociSTREAMitfsi
! if(mod(TStep,100).eq.0)then
    casename = 'TecDef'
    outputFileNLAMS=casename(1:6)//achar(int(TStep/100000)+48)//achar(int(mod(TStep,100000)/10000)+48)//&
&achar(int(mod(TStep,10000)/1000)+48)//achar(int(mod(TStep,1000)/100)+48)//achar(int(mod(TStep,100)/10)+48)&
&//achar(mod(TStep,10)+48)//'.dat'
    unitnum = 1000
       open(unit=unitnum,file=outputFileNLAMS,action='write')
       rewind(unitnum)
          write(unitnum,*)'title = "Sample finite-element data"'
          write(unitnum,*)'VARIABLES = "X", "Y", "Z"'
          write(unitnum,*)'ZONE F=FEPOINT, ET=TRIANGLE, N=',NumNodes,',E=',NumElems
       do i=1,NumNodes
          write(unitnum,'(3(f18.9,1x))')CurrCoordInert(i,1:3)
       enddo
        do i=1,NumElems
           write(unitnum,'(i6,1x,i6,1x,i6)') Connect(i,1:3)
        enddo
       close(unitnum)
! endif
endif
! if (LociSTREAMRank==0)  print *,   'A5', LociSTREAMitfsi
! !End of the newly added stuff

 return 
end subroutine nlams_output

subroutine nlams_locistream_output_step (T_step,Totdisp)
use mesh_connect_vars, only: NumDofPerNode,NumNodes,Connect,UndefGrid,NumElems
!use nlams_init,         only: AnsSize !commented out for including tecplot
!files output at each step
use nlams_init,         only: AnsSize,DefGrid
use nlams_undeformed,   only: T,MemOri
use shared_time,        only:TheTime

!Local variables
integer :: kstart,kend,i,unitnum
real(8),dimension(3,1)::Aux,Tdisp1
character(len=80)::outputFileNLAMS,casename

!intent(in) variables
real(8),dimension(AnsSize,1),intent(in)::TotDisp
integer,intent(in)                     ::T_step

write(393,*)MemOri(1:3,1)
!721,722,723,121 are the numbers for the old mesh
   !Write ASCII output file for tip displacement at the current load level
     Aux(1,1) = Totdisp(7,1)
     Aux(2,1) = Totdisp(8,1)
     Aux(3,1) = Totdisp(9,1)
     Aux(1:3,1) = Aux(1:3,1) + UndefGrid(2,1:3) + MemOri(1:3,1)

     Tdisp1(:,1) = matmul(T,Aux(:,1))

     write(204,*)TheTime,Tdisp1(3,1)

!Newly added to generate tecplot files at each time-step
    casename = 'TecDef'
    outputFileNLAMS=casename(1:6)//achar(int(T_step/1000)+48)//achar(int(mod(T_step,1000)/100)+48)&
&//achar(int(mod(T_step,100)/10)+48)//achar(mod(T_step,10)+48)//'.dat'

    unitnum = 1000

    open(unit=unitnum,file=outputFileNLAMS,action='write')
    !Write TECPLOT unstructured FE data file out
          write(unitnum,'(A)')'title = "Sample finite-element data"'
          write(unitnum,'(A)')'variables = "x", "y", "z"'
          write(unitnum,'(A,I5,A,A,I5,A,A)')'zone n=',NumNodes,',','e=',NumElems,',','DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
       do i=1,NumNodes
           write(unitnum,*)DefGrid(i,1:3)
       enddo

        do i=1,NumElems
           write(unitnum,*)Connect(i,1:3)
        enddo
        close(unitnum)
!End of the newly added stuff

  return
 end subroutine nlams_locistream_output_step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module nlams_interf
