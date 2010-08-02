!-> Copyright by The University of Michigan, Aerospace Department. 2004
!
!-> Module.- NLAMS_HEADER. Rafa Palacios. 9Sep2002 / 29Sep2004
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!    Open general I/O Files and write header on them.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_header
 use mod_shared
 implicit none
!
! (No Public variables in this module).
!
 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine nlams_header_main
 use LociSTREAMdata, only: LociSTREAMRank
 
 implicit none
!
! Internal flag to activate structural analysis modules.
!

! Output file names
  if (LociSTREAMRank==0)  then
  open(204 ,file='TipDispLoadLevel.dat',action='write')
!  open(339,file='Energy.dat',action='write')
  open(2302,file='CSDDelta.dat',action='write')
  open(393 ,file='MemOri.dat',action='write')
  open(340 ,file='PrescribedMotionFlap.dat',action='write')
  open(345 ,file='PrescribedMotionPlunge.dat',action='write')
  endif
! Not using
!  open(2403,file='FinalAnsCSD.dat',action='write')
!  open(202,file='output.dat',action='write')

  StructAnalysis=.true.
!
! Open the file with inputs in MSC/NASTRAN formatting.
!
   InputFile = 'nlamstest'
!   open (unit= iuInput,                &
!&      file= trim(InputFile)//'.nlam', &
!&      action='READ')
!
! Open main output file.
!
!  open (unit= iuOutput,                &
!&      file= trim(InputFile)//'.nout', &
!&      action= 'WRITE')
!
! Write header on the screen.
!
if (LociSTREAMRank==0)  print *, 'Starting NLAMS, Release 1.1'
!
! Write header in the General Output file.
!
!  write (iuOutput,'(9(40X,A,/))')                           &
!&      '*************************************************', &
!&      '*                                               *', &
!&      '*            The University of Michigan         *', &
!&      '*      Department of Aerospace Engineering      *', &
!&      '*                                               *', &
!&      '*          NonLineAr Membrane Shell Solver      *', &
!&      '*                Release 1.1                    *', &
!&      '*                                               *', &
!&      '*************************************************'

  return
 end subroutine nlams_header_main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module nlams_header

