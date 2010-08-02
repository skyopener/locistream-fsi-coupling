!-> Copyright by The University of Michigan, Aerospace Department. 2003
!
!-> Module.- MOD_SHARED
!
!-> Description.-
!
!  This module defines the set of parameters and the error subroutine
!  that are shared by the subprograms in program SmartRotor.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mod_shared
 implicit none
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! Parameters for Preallocated memory (only for aero part).
!
  integer,parameter:: MaxBlades= 2  ! Maximum Number of Blades.
  integer,parameter:: MaxElems= 20  ! Maximum Number of Elements.
!
! Logical Units.
!
  integer,parameter:: iuAllOut= 0    ! Send data to all output files.
  integer,parameter:: iuInpScr= 5    ! Default input file (screen).
  integer,save     :: iuOutScr       ! Default output file (screen).
  integer,parameter:: iuInput =11    ! Input File  (.nlam)
  integer,parameter:: iuOutput=12    ! Output file (.nout)
  integer,parameter:: iuDat3  =13    ! Tecplot output file (_grid.dat)
  integer,parameter:: iuDat6  =14    ! Tecplot output file (_aero.dat)
  integer,parameter:: iuN3g   =15    ! Output file (.n3g)
  integer,parameter:: iuN3d   =16    ! Output file (.n3d)
  integer,parameter:: iuN3s   =17    ! Output file (.n3s)
  integer,parameter:: iuTmp1  =31    ! Temporary file.
!
! Real Constants.
!
  integer,parameter:: MaxNestedCords=99              ! Maximum nested coordinate definitions.
  real(8),parameter:: MaxTime=9.d99                  ! Maximum time allowed in the code.
  real(8),parameter:: EpsilonDefault=1.d-8           ! Small number
  real(8),parameter:: Pi =3.14159265358979323846264338327950288 ! Pi
  real(8),parameter,dimension(3,3):: Unit= &         ! Unit matrix.
&         reshape((/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/),(/3,3/))
!
! Error codes.
!
  integer,parameter:: ErrNeededLab    =   1
  integer,parameter:: ErrRepeatLab    =   2
  integer,parameter:: ErrInvaldLab    =   3
  integer,parameter:: ErrWrongData    =   4
  integer,parameter:: ErrBadFile      =   5
  integer,parameter:: ErrAllocMem     =   6
  integer,parameter:: ErrSolver       =   7
  integer,parameter:: ErrInputFile    =   9
  integer,parameter:: ErrBadMesh      =  10
  integer,parameter:: ErrConverg      = 111
  integer,parameter:: ErrInvalidParam =  12
  integer,parameter:: ErrNLConverg    =  13
  integer,parameter:: ErrMdiceAlloc   =  14
  integer,parameter:: ErrMdiceGeneral =  15
  integer,parameter:: ErrNotSteady    =  16
  integer,parameter:: ErrSmallMemory  =  17
  integer,parameter:: ErrEigSol       =  18
!
! Input file and flags for modules.
!
 character(len=256),save:: InputFile  ! Name of the input file.
 logical,save::        StructAnalysis ! True if structural module is included.
 logical,save::        AerodAnalysis  ! True if aero module is included.
 logical,save::        MdiceAnalysis  ! True if running into MDICE.
!
! Buffer for reading input file lines.
!
 character(len=10000),save:: DataBlock     ! Data read in an entry.

 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ERROR.
!
!->Description.-
! 
!  Handle Error in Program Execution.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine error (iOut, CodErr, TxtErr, ID)
  implicit none
!
!-> Input/Output Variables.
!
  integer, intent(in) :: iOut    ! Logic Unit of the output file.
  integer, intent(in) :: CodErr  ! Error codes as defined in params of module.
  character(len=*), intent(in) :: TxtErr  ! Text to be written.
  integer, intent(in), optional:: ID      ! ID of an element with error.
!
!-> Local variables.
!
  integer              :: i        ! Counter.
  integer, allocatable :: Files(:) ! Logic Unit where data will be written.
!
! Set the files where the text error is going to be written.
!
  if (iOut .eq. iuAllOut) then
    allocate(Files(2))
    Files(1)=iuOutScr
    Files(2)=iuOutput
  else
    allocate(Files(1))
    Files(1)=iOut
  end if
!
! Write the information on the selected output file.
!
  do i=1, size(Files)
    write (Files(i),'(/A5,$)') '---> '
!
    select case (CodErr)



    case (ErrNeededLab)
     write (Files(i),'(A,$)') 'ERROR: Label ' // adjustl(TxtErr) // &
&                             ' was not found in Input File.'
     if (present(ID)) then
       write (Files(i),'(A,I8)') ' with ID', ID
     else
       write (Files(i),'(A)') ' '
     end if



    case (ErrRepeatLab)
     write (Files(i),'(A)') 'ERROR: Label ' // adjustl(TxtErr) // &
&                           ' cannot be repeated in Input File.'


    case (ErrInvaldLab)
     write (Files(i),'(A,$)') 'ERROR: Invalid Format in Label '// &
&                             adjustl(TxtErr)
     if (present(ID)) then
       write (Files(i),'(A,I8)') ' with ID', ID
     else
       write (Files(i),'(A)') ' '
     end if



    case (ErrWrongData)
     write (Files(i),'(A)') 'ERROR: Incorrect data in the Input File:'
     write (Files(i),'(12X,A,$)') adjustl(TxtErr)
     if (present(ID)) then
       write (Files(i),'(A,I8)') ' with ID', ID
     else
       write (Files(i),'(A)') ' '
     end if



    case (ErrBadFile)
     write (Files(i),'(A)') 'ERROR: Error on File '// adjustl(TxtErr)



    case (ErrAllocMem)
     write (Files(i),'(A)') 'ERROR: Not enough memory to allocate '// &
&                           'variable '//adjustl(TxtErr)



    case (ErrSolver)
     write (Files(i),'(A)') 'ERROR: Error in the linear equations solver.'
     write (Files(i),'(12X,A)') adjustl(TxtErr)



    case (ErrInputFile)
     write (iuOutScr,'(A)') 'ERROR: Input file not found: '//adjustl(TxtErr)



    case (ErrBadMesh)
     write (Files(i),'(A)') 'ERROR: Problem found in mesh definition.'
     write (Files(i),'(14X,A,$)') adjustl(TxtErr)
     if (present(ID)) then
       write (Files(i),'(A,I8)') ' with ID', ID
     else
       write (Files(i),'(A)') ' '
     end if


    case (ErrConverg)
     write (Files(i),'(A)') 'WARNING: Convergence Error'
     write (Files(i),'(14X,A)') adjustl(TxtErr)


    case (ErrInvalidParam)
     write (Files(i),'(A)') 'ERROR: Invalid Format in Entry nPARAM '// &
&                             adjustl(TxtErr)


    case (ErrNLConverg)
     write (Files(i),'(A)') 'ERROR: Max Number of iterations was reached'
     write (Files(i),'(14X,A)') adjustl(TxtErr)



    case (ErrMDICEAlloc)
     write (Files(i),'(A)') 'ERROR: MDICE could not allocate memory.'
     write (Files(i),'(14X,A)') adjustl(TxtErr)



    case (ErrMDICEGeneral)
     write (Files(i),'(A)') 'ERROR: The MDICE controller has returned an error.'
     write (Files(i),'(14X,A)') adjustl(TxtErr)


    case (ErrNotSteady)
     write (Files(i),'(A)') 'ERROR: No steady-state solution is possible.'
     write (Files(i),'(14X,A)') adjustl(TxtErr)


    case (ErrSmallMemory)
     write (Files(i),'(A,$)') 'ERROR: Not enough memory for matrix:'
     write (Files(i),'(2X,A)') adjustl(TxtErr)


    case (ErrEigSol)
     write (Files(i),'(A)') 'ERROR: Error in the eigenvalue solver.'
     write (Files(i),'(14X,A)') adjustl(TxtErr)

    end select
  end do
!
! Codes lower than 100 implies Stop of program execution.
!
  if (CodErr .lt. 100) then
   write (*,'(/A)') ' Program NLAMS terminates abnormally.'
!   pause 'Press RETURN to continue'
   stop ' '
  end if
!
  deallocate (Files)
!
  return
 end subroutine error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module mod_shared

