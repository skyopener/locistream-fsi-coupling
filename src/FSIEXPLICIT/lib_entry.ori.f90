!
!->Module LIB_ENTRY. Rafa Palacios. 23Ago2002
!
!->Description.-
!
! Read/Write Data lines identified by an entry code in the input file.
! Fixed format of 8 characters by field is used.
!
!->Module Subroutines.-
!
! check_wrong_entry: Check wrong entries in the input file.
! get_next_entry:    Get next occurrence of a given entry.
! get_number_entry:  Get number of occurrences of a given entry.
! write_entry:       Write a data line in fixed format.
!
!->Internal Subroutines:
!
! get_entry:           Read a card in fixed format.
! store_entry_in_list: Store read entry labels in a DataBase.
! write_real_in_entry_field: Write real output with accuracy.
!
!->Remarks:
!
! 1) This module should be compiled with Case-Sensitive options 'on'.
!
! 2) External subroutines used in the module:
!
!     * str_upcase: Convert a string chain to Capitals.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_entry
 implicit none
!
!-> Private variables to the library:
!   --------------------------------
!
! List of labels that were found in the input file.
 character (len=8), save, private, allocatable, dimension(:) :: LisLabels
!
! Number of labels in the list.
 integer, save, private :: NumLabels=0
! 
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine GET_NUMBER_ENTRY
!
!-> Description.- 
! 
!  Get the number of occurrences of a given Entry label in a file.
!  The file must come in MSC/NASTRAN fixed input format and the Entry
!  label is defined in the first 8 characteres of the Entry.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine get_number_entry (iuInp,Label,NumEntries,MaxFields)
  use lib_tools
!
!-> Input variables.
!
  integer,         intent(in) :: iuInp      ! Logical iuInp for input.
  character(len=8),intent(in) :: Label      ! Identifying label.
  integer,         intent(out):: NumEntries ! Total # of occurrences of Label.
  integer,optional,intent(out):: MaxFields  ! Max number of data fields.
!
!-> Local variables.
!
  character (len=6400) :: DataBlock  ! Dummy output text.
  character (len=8)    :: LabelRead  ! Read value for the label.
  integer :: iError                  ! Error code. its value is non-zero when
                                     ! a label was not found.
  integer :: NumFields               ! Number of fields in an entry.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Go to the first line of the input file and set initial values of variables.
!
  rewind (iuInp)
  iError=0
  NumEntries=0
  if (present(MaxFields)) MaxFields=0
!
! Read new entries until the given label is found.
!
  do while (iError.eq.0)
    call get_entry (iuInp,LabelRead,DataBlock,iError)

    if ( (iError.eq.0) .and. (LabelRead.ne.' ') ) then
      call str_upcase (LabelRead, LabelRead)
      if (adjustl(LabelRead).eq.Label) then

        NumEntries= NumEntries+1
!
! Get the number of DataBlock fields if NumFields was required.
!
        if (present(MaxFields)) then
          if (index(Label,'*').eq.len_trim(Label)) then
            NumFields= ceiling(real(len_trim(DataBlock))/16.)
          else
            NumFields= ceiling(real(len_trim(DataBlock))/8.)
          end if
          if (NumFields.gt.MaxFields) MaxFields=NumFields
        end if

      end if
    end if
  end do
!
! Store the label in the label list and the output value.     
!
  if (NumEntries .gt. 0) then
    call store_entry_in_list (label)
  end if
!
! End of subroutine.
!
  rewind (iuInp)
  return
 end subroutine get_number_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine.- GET_NEXT_ENTRY
! 
!->Description.- 
! 
!  Read information for the next aparition of a given Entry in a file.
!  The file must come in MSC/NASTRAN fixed input format and the Entry code
!  is defined in the first 8 characteres of the Entry.
!
!-> Remarks.-
!
!  1) Variable DATA includes only Fields 2 to 9 (2-5 in long field format)
!     in all input lines of the current Entry. The size of the variable 
!     must be determined by the calling subroutine, and since all data 
!     fields in a line are stored (8x8 or 4x16 fields per line), the size
!     of DATA in the calling subroutine must be a multiple of 64.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine get_next_entry (iuInp, Label, DataBlock, iError, NumFields)
  use lib_tools
!
!-> Input/Output Variables.
!
  integer,           intent(in)  :: iuInp      ! Logical unit for input.
  character (len=8), intent(in)  :: Label      ! Identifying Label.
  character (len=*), intent(out) :: DataBlock  ! Output text.
  integer,           intent(out) :: iError     ! Error code.
  integer, optional, intent(out) :: NumFields  ! Number of data fields.
!
!-> Local Variable.
!
  logical :: LabelFound          ! Flag.
  character (len=8) :: LabelInp  ! Read value for the label.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read new entries until the given label is found.
!
  DataBlock=' '
  LabelFound= .false.
  do while (.not. LabelFound)
    call get_entry (iuInp, LabelInp, DataBlock, iError)
!
! Convert the label code to capitals and compare with the given one.
!
    if (len_trim(LabelInp) .gt. 0 ) then       ! Black labels not allowed
      call str_upcase (LabelInp, LabelInp)
      if (adjustl(LabelInp) .eq. Label) LabelFound=.true.
    end if
!    
  end do
!
! Get the number of DataBlock fields if NumFields was required.
!
  if (present(NumFields)) then
    if (index(Label,'*') .eq. len_trim(Label)) then
      NumFields= ceiling(real(len_trim(DataBlock))/16.)
      else
      NumFields= ceiling(real(len_trim(DataBlock))/8.)
    end if
  end if
!
  return
 end subroutine get_next_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine CHECK_WRONG_ENTRY.
! 
!->Description.- 
! 
!  Check if the input file in MSC/NASTRAN format has entries different
!  to the set store in the Entries List.
!  
!  The Entries List was generated each time the main program called the
!  subroutine get_number_entries.
!
!-> Remarks.-
!
!  1) Outputs are warnings on the Output File/Screen with the non-used lines.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine check_wrong_entry (iuInp, iuOut)
  use lib_tools
!
!--> Input variables.
!
  integer, intent(in) :: iuInp   ! Logical unit of input file.
  integer, intent(in) :: iuOut   ! Logical unit of output file/Screen.
!
!--> Internal variables.
!
  character (len=1) :: DataBlock ! Dummy entry data block (not stored).
  logical :: LabelFound          ! Flag for a label in the list.
  integer :: i                   ! Counter.
  integer :: iError              ! Read error code.
  integer :: NumWrongLabels      ! Num of wrong labels in the input file.
  character (len=8) :: Label     ! Label of the entry.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Initialization of variables.
!
  rewind (iuInp)
  iError=0
  NumWrongLabels=0
!
! Read all entries in the input file and check if the label of each entry
! is included in the label list.
!
  do while (iError .eq. 0)
    LabelFound=.false.
    call get_entry (iuInp, Label, DataBlock, iError)
!
! Blank lines are not given back by Get_Entry, therefore, if a line has
! blank label it must have non-blank data block.
!
    if (len_trim(Label) .gt. 0) then
      call str_upcase (Label, Label)
      do i=1, NumLabels
        if (LisLabels(i) .eq. label) LabelFound=.true.
      end do
    end if
!
! Write the error code on screen when a non defined label was found.
!
    if (.not.LabelFound) then
      NumWrongLabels=NumWrongLabels+1
!      
      if (NumWrongLabels.eq.1) then
        write (iuOut,'(/a/)') &
&         '***************************************************************'
        write (iuOut,'(a,/)') &
&         '     Warning!! These input cards were not recognised:'
      end if
!
      if (len_trim(Label) .gt. 0) then
        write (iuOut,'(25x,a)') Label
      else
        write (iuOut,'(26x,a)') '(entry without label)'
      end if
!
    end if
  end do
!
  if (NumWrongLabels .ge. 1) then
    write (iuOut,'(/a,i3/)') &
&      '     number of wrong labels: ',NumWrongLabels
    write (iuOut,'(a)') &
&      '***************************************************************'
  end if
!  
  return
 end subroutine check_wrong_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine.- WRITE_ENTRY
! 
!->Description.- 
! 
!  Write fields if MSC.Nastran fixed format.
!
!-> Remarks.-
!
!  1) Short fields (8-characters length) is used.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine write_entry (iuOut, Label, i1,i2,i3,i4,i5,i6,i7,i8, &
&                                      r1,r2,r3,r4,r5,r6,r7,r8, &
&                                      d1,d2,d3,d4,d5,d6,d7,d8, &
&                                      c1,c2,c3,c4,c5,c6,c7,c8, Cont)
!
!-> I/O Variables.
!
  integer               , intent(in) :: iuOut
  character(*), optional, intent(in) :: Label
  integer,      optional, intent(in) :: i1,i2,i3,i4,i5,i6,i7,i8
  real(4),      optional, intent(in) :: r1,r2,r3,r4,r5,r6,r7,r8
  real(8),      optional, intent(in) :: d1,d2,d3,d4,d5,d6,d7,d8
  character(*), optional, intent(in) :: c1,c2,c3,c4,c5,c6,c7,c8
  logical,      optional, intent(in) :: Cont
!
!-> Local Variables.
!
  character(len=80) :: LineOut    ! Line to be written in the output file.
!
  LineOut=' '
!
! Write the label.
!
  if (present(Label)) then
    LineOut(1:8)=adjustl(Label)
  else
    LineOut(1:1)='+'
  end if
!
! Write the 1st field.
!
  if (present(i1)) then
    write (LineOut(9:16),'(I8)') i1
  else if (present(r1)) then
    call write_real_in_entry_field (LineOut(9:16), dble(r1))
  else if (present(d1)) then
    call write_real_in_entry_field (LineOut(9:16), d1)
  else if (present(c1)) then
    LineOut(9:16)= adjustr(c1)
  end if
!
! Write the 2nd field.
!
  if (present(i2)) then
    write (LineOut(17:24),'(I8)') i2
  else if (present(r2)) then
    call write_real_in_entry_field (LineOut(17:24), dble(r2))
  else if (present(d2)) then
    call write_real_in_entry_field (LineOut(17:24), d2)
  else if (present(c2)) then
    LineOut(17:24)= adjustr(c2)
  end if
!
! Write the 3rd field.
!
  if (present(i3)) then
    write (LineOut(25:32),'(I8)') i3
  else if (present(r3)) then
    call write_real_in_entry_field (LineOut(25:32), dble(r3))
  else if (present(d3)) then
    call write_real_in_entry_field (LineOut(25:32), d3)
  else if (present(c3)) then
    LineOut(25:32)= adjustr(c3)
  end if
!
! Write the 4th field.
!
  if (present(i4)) then
    write (LineOut(33:40),'(I8)') i4
  else if (present(r4)) then
    call write_real_in_entry_field (LineOut(33:40), dble(r4))
  else if (present(d4)) then
    call write_real_in_entry_field (LineOut(33:40), d4)
  else if (present(c4)) then
    LineOut(33:40)= adjustr(c4)
  end if
!
! Write the 5th field.
!
  if (present(i5)) then
    write (LineOut(41:48),'(I8)') i5
  else if (present(r5)) then
    call write_real_in_entry_field (LineOut(41:48), dble(r5))
  else if (present(d5)) then
    call write_real_in_entry_field (LineOut(41:48), d5)
  else if (present(c5)) then
    LineOut(41:48)= adjustr(c5)
  end if
!
! Write the 6th field.
!
  if (present(i6)) then
    write (LineOut(49:56),'(I8)') i6
  else if (present(r6)) then
    call write_real_in_entry_field (LineOut(49:56), dble(r6))
  else if (present(d6)) then
    call write_real_in_entry_field (LineOut(49:56), d6)
  else if (present(c6)) then
    LineOut(49:56)= adjustr(c6)
  end if
!
! Write the 7th field.
!
  if (present(i7)) then
    write (LineOut(57:64),'(I8)') i7
  else if (present(r7)) then
    call write_real_in_entry_field (LineOut(57:64), dble(r7))
  else if (present(d7)) then
    call write_real_in_entry_field (LineOut(57:64), d7)
  else if (present(c7)) then
    LineOut(57:64)= adjustr(c7)
  end if
!
! Write the 8th field.
!
  if (present(i8)) then
    write (LineOut(65:72),'(I8)') i8
  else if (present(r8)) then
    call write_real_in_entry_field (LineOut(65:72), dble(r8))
  else if (present(d8)) then
    call write_real_in_entry_field (LineOut(65:72), d8)
  else if (present(c8)) then
    LineOut(65:72)= adjustr(c8)
  end if
!
! Write the continuation mark.
!
  if (present(Cont)) then
    if (Cont) LineOut(80:80)='+'
  end if
!
! Copy to the output file.
!
  write (iuOut,'(A80)') LineOut
!
  return
 end subroutine write_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function.- STORE_ENTRY_IN_LIST
!
!-> Description.-
!
! Store the identification label of an entry in the list of read values.
!
!-> Input Variables.-
!
! Label : Label in the selected Entry.
!
!-> Output of the function.-
!
!   The number of ocurrences of the entry.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine store_entry_in_list (Label)
!
!-> Input variables.
!
  character (len=8), intent(in) :: Label
!
!-> Local variables.
!
! Variable used to temporarily store the label list.
  character (len=8), allocatable, dimension(:) :: LisAux
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First label in the list.
!
  if (NumLabels .eq. 0) then
    allocate (LisLabels(1))
    LisLabels(1)=Label
!
! Subsequent labels are added at the end of the list. Since the dimension
! of the list must increase, it is first stored in the auxiliar array.
!
  else
    allocate (LisAux(NumLabels))
    LisAux=LisLabels
    deallocate (LisLabels)
    allocate (LisLabels(NumLabels+1))
    LisLabels(1:NumLabels)=LisAux
    LisLabels(NumLabels+1)=Label
    deallocate (LisAux)
  end if
  NumLabels=NumLabels+1
!
  return
 end subroutine store_entry_in_list
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!->Subroutine.- GET_ENTRY
! 
!->Description.- 
! 
!  Reads all fields of a new entry in MSC/NASTRAN Format. Fields have a length
!  of 8 characters.
!
!-> Remarks.-
!
!  1) Data for an Entry can be given in several lines by adding symbol '+' at 
!     Field #10 for small field format and '*' for long field format.
!
!  2) Field #1 in the first line is the identifying label. In continuation
!     lines this field is skipped.
!
!  3) Long field format entries are identified by a label ending with '*'.
!
!  4) Commented lines begin with '$'.
!
!  5) If iuInp==0, data is read from screen.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine get_entry (iuInp, Label, DataBlock, iError)
!
!-> Input/Output variables.
!
  integer,           intent(in)  :: iuInp      ! Logic unit for input file.
  character (len=8), intent(out) :: Label      ! Identifying Label.
  character (len=*), intent(out) :: DataBlock  ! Fields 2-9 of input line.
  integer,           intent(out) :: iError     ! Error code.
!
!-> Local variables.
!
  character(len= 1) :: ContSymbol ! Continuation symbol.
  character(len=80) :: LineInp    ! Input data line.
  integer :: NumLines             ! Number of lines already read.
  logical :: ReadNewLine          ! Flag to read (t) or not (f) a new line.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Initialization of variables.
!
  NumLines=0
  ReadNewLine=.true.
  do while (ReadNewLine)
!
! Read the next line.
!
    LineInp=' '
    if (iuInp.eq.0) then
      read (*,fmt='(a80)',iostat=iError) LineInp
    else
      read (iuInp,fmt='(a80)',iostat=iError) LineInp
    end if
    if (iError .ne. 0) return     ! Exit if an error was found.
!
! Check if the new line is commented.
!
    if ((LineInp(1:1).eq.'$') .or. (len_trim(LineInp).eq.0)) then
      ReadNewLine=.true.
    else
!
! For the first line, the label is got from the first field.
!
      NumLines=NumLines+1
      if (NumLines .eq. 1) then
        Label= adjustl(LineInp(1:8))
        if (index(Label,'*') .eq. len_trim(Label)) then
          ContSymbol='*'
        else
          ContSymbol='+'
        end if
      end if
!
! Get the DataBlock. This subroutine stores data up to the maximum length
! given by the user and does not give an error code if data are out of these
! limits.
!
      if (NumLines*64 .le. len(DataBlock)) then
        DataBlock ((NumLines-1)*64+1:)= LineInp(9:72)
      end if
!
! Handle continuation symbol at field #10.
!
      if (index(adjustl(LineInp(73:80)),ContSymbol) .eq. 1) then
        ReadNewLine=.true.
      else
        ReadNewLine=.false.
      end if
!
    end if
  end do
!
  return
 end subroutine get_entry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!->Subroutine.- WRITE_REAL_IN_ENTRY_FIELD
! 
!->Description.- 
! 
!  Write a REAL*8 number in a 8-characters field using the optimum format
!  for accuracy.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine write_real_in_entry_field (CharOut, RealIn)
!
!-> I/O Variables.
!
  character(len=8), intent(out) :: CharOut
  real(8),          intent(in)  :: RealIn
!
! Local Variables.
!
  real(8) :: ExpNum   ! Integer Logarithm of RealIn.
!
  CharOut=' '
  if (abs(RealIn).le.1.d-10) then
    CharOut(7:8)='0.'
  else
    ExpNum=dint(dlog10(dabs(RealIn)))
!
    if (ExpNum.ge.100.) then
	  if (RealIn.ge.0.) then
        CharOut(2:8)= '9.99+99'
      else
        CharOut(1:8)='-9.99+99'
      end if
!
    else if (ExpNum.ge.10. .and. ExpNum.lt.100.) then
      if (RealIn.ge.0.) then
        write (CharOut(1:5),'(F5.3)') RealIn*(1.d1**(-ExpNum))
      else
        write (CharOut(1:5),'(F5.2)') RealIn*(1.d1**(-ExpNum))
      end if
      CharOut(6:6)='+'
      write (CharOut(7:8),'(I2)') int(ExpNum)
!
    else if (ExpNum.ge.0. .and. ExpNum.lt.10.) then
      if (RealIn.ge.0.) then
        write (CharOut(1:6),'(F6.4)') RealIn*(1.d1**(-ExpNum))
      else
        write (CharOut(1:6),'(F6.3)') RealIn*(1.d1**(-ExpNum))
      end if
      CharOut(7:7)='+'
      write (CharOut(8:8),'(I1)') int(ExpNum)
!
    else if (ExpNum.gt.-10. .and. ExpNum.lt.0.) then
      if (RealIn.ge.0.) then
        write (CharOut(1:6),'(F6.4)') RealIn*(1.d1**(-ExpNum))
      else
        write (CharOut(1:6),'(F6.3)') RealIn*(1.d1**(-ExpNum))
      end if
      write (CharOut(7:8),'(I2)') int(ExpNum)
!
    else if (ExpNum.gt.-100. .and. ExpNum.le.-10.) then
      if (RealIn.ge.0.) then
        write (CharOut(1:5),'(F5.3)') RealIn*(1.d1**(-ExpNum))
      else
        write (CharOut(1:5),'(F5.2)') RealIn*(1.d1**(-ExpNum))
      end if
      write (CharOut(6:8),'(I3)') int(ExpNum)
    end if
  end if
!
 return
 end subroutine write_real_in_entry_field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_entry
