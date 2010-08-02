!->Copyright by The University of Michigan, Aerospace Department. 2002
!
!->Module LIB_TOOLS. Rafa Palacios. 11Aug2002
!
!->Description.-
!
!  This module includes different elementary tools.
!
!->Subroutines:
!
!   add_blockmat:  Add block matrix to bigger matrix.
!   add_blockvec:  Add "block" vector to bigger vector.
!   id_extract:    Extract order of ID in a List.
!   new_line:      Write a blank line in a file.
!   str_upcase:    Convert a string to capitals.
!   under_line:    Echo of computations in output file.
!   write_comment: Write commented line in output file.
!   write_header:  Write header in columns.
!   write_title:   Write main title.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_tools
 implicit none
!
!   (There are no public variables in this module).
!
 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ADD_BLOCKMAT.
!
!->Description.-
!
!   Add block matrix to a bigger matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine add_blockmat(Aglobal,i0,j0,Alocal)

  real(8),intent(inout):: Aglobal(:,:)   ! Global matrix
  integer,intent(in)   :: i0             ! Initial row.
  integer,intent(in)   :: j0             ! Initial column.
  real(8),intent(in)   :: Alocal(:,:)    ! Local matrix.

  integer :: M,N   ! Dimension of Alocal

  M=size(Alocal,DIM=1)
  N=size(Alocal,DIM=2)

  Aglobal(i0+1:i0+M,j0+1:j0+N)= Aglobal(i0+1:i0+M,j0+1:j0+N) + Alocal

  return
 end subroutine add_blockmat



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ADD_BLOCKVEC.
!
!->Description.-
!
!   Add block matrix to a bigger matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine add_blockvec(Aglobal,i0,Alocal)

  real(8),intent(inout):: Aglobal(:)     ! Global vector
  real(8),intent(in)   :: Alocal(:)      ! Local vector.
  integer,intent(in)   :: i0             ! Initial row.

  integer :: M   ! Dimension of Alocal
  M=size(Alocal,DIM=1)

  Aglobal(i0+1:i0+M)= Aglobal(i0+1:i0+M) + Alocal

  return
 end subroutine add_blockvec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Function FACTORIAL.
!
!->Description.-
!
!   Factorial of a integer.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real(8) function dfactorial (n)

  integer,intent(in) :: n
!
  integer:: k,fact

  if (n.eq.0) then
    fact=1

  else
    k=n
    fact=k
    do while (k.gt.1)
      k=k-1
      fact=fact*k
    end do
  end if

  dfactorial=1.d0*fact
  return
 end function dfactorial


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine ID_EXTRACT
!
!->Description.-
!
!   Search the label ID_Input in the ID_List and returns the index of
!   the first appearance of the label in ID_List (zero if not present).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 integer function id_extract (ID_Input,ID_List)
!
!-> I/O variables.
!
  integer,intent(in):: ID_Input    ! Input ID label.
  integer,intent(in):: ID_List(:)  ! List where it will be searched.
!
! Local Variables.
!
  integer :: NumEl ! Number of elements in the list.
  integer :: i     ! Counter.
!
! Extract the position of the given label in the DataBase.
!
  id_extract=0
  do i=1,size(ID_List)
    if (ID_Input.eq.ID_List(i)) then
      id_extract=i
      exit
    end if
  end do
!
  return
 end function id_extract





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine NEW_LINE.
!
!->Description.-
!
!   Output a new line to a file connected to unit out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine new_line(out)
!
  integer,intent(in)::out
  write(out,'(1x,/)')
  return
 end subroutine new_line






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine STR_UPCASE.
!
!->Description.-
!
! Convert to Capitals a given string.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine str_upcase (StrUpp,StrInp)
!
!-> Input/Output Variables.
!
  character(len=*), intent(in) :: StrInp     ! Input String including lower cases.
  character(len=*), intent(out):: StrUpp     ! Output String in Capitlals.
!
!-> Local Variables.
!
  integer, parameter :: TotChar=26           ! Number of allowed characters.
  integer :: i,j                             ! Counters.
!
  character (len=1), dimension(TotChar) :: AdmUpp ! Allowed characts in Capts.
  character (len=1), dimension(TotChar) :: AdmLow ! Allowed characts in lower.
!
  data AdmUpp/'A','B','C','D','E','F','G','H','I','J','K','L','M',  &
&              'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/
  data AdmLow/'a','b','c','d','e','f','g','h','i','j','k','l','m',  &
&              'n','o','p','q','r','s','t','u','v','w','x','y','z'/
!
! Substitution in the input string.
!
  StrUpp(1:len_trim(StrInp))= trim(StrInp)
  do i= 1, len_trim(StrInp)
    do j=1, TotChar
      if (StrInp(i:i) .eq. AdmLow(j)) StrUpp(i:i)=AdmUpp(j)
    end do
  end do
!
  return
 end subroutine str_upcase






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine UNDER_LINE.
!
!->Description.-
!
!   Output a line with a bunch of '=' to  out 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine under_line(out)
!  
  integer,intent(in)::out
  write(out,'(1x,50("="))')
!  
  return
 end subroutine under_line






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!->Subroutine WRITE_COMMENT.
!
!->Description.-
!
!   Write a commented line in the given output file.
!   Commented lines start with the $ symbol.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine write_comment (iuOut,Text,Line,Date)
!
!-> I/O Variables.
!
  integer, intent(in) ::           iuOut    ! Logical unit of the file.
  character(len=*), intent(in) ::  Text     ! Text to be written.
  logical, intent(in), optional :: Line     ! Option to write a line of '$'.
  logical, intent(in), optional :: Date     ! Option to write the date.
!
  write (iuOut,'(A)') '$'
  write (iuOut,'(A)') '$ '//Text
  write (iuOut,'(A)') '$'
  if (present(Line)) then
    if (Line) write (iuOut,'(A80)') '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'// &
&                                   '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  end if
!
  return
 end subroutine write_comment
 
 

 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine WRITE_HEADER
!
!-> Description.-
!
!     Writes headers for the columns in the output file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine write_header (File, NumCols, Headers)
!
!-> Input Variables.
!
  integer      File        ! Output file.
  integer      NumCols     ! Number of columns.
  character*8  Headers(*)  ! Headers for the columns.
!
!-> Local variables.
!
  integer      i
!
  write (File, '(/16(4X,A8,4X))') (Headers(i), i=1, NumCols)
!
  return
 end subroutine write_header





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine WRITE_TITLE
!
!-> Description.-
!
!     Writes a main title in a given file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine write_title (File, Title)
!
!-> Input Variables.
!
  integer      File        ! Output file.
  character(*) Title       ! Title to be written.
!
!-> Local variables.
!
  integer      i
!
!
  write (File,'(//10X,A)') Title
  write (File,'(10X,100A1)') ('=',i=1,len_trim(Title)+1)
!
  return
 end subroutine write_title

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lib_tools
 
