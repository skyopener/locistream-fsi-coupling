!-> Subroutine SHARED_HEADER. Rafa Palacios. 29Mar2003
!   Modified: Satish Chimakurthi Sep 2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!    Get the generic name of the input file from screen and
!    write headers.
!
!-> Remarks.-
!
!  1) The use of the non-standard subroutine GETARG with the Digital 
!     Compiler requires the reference to the module dflib.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shared_header
 use mod_shared
 implicit none

 contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine shared_header_main
  !use dflib
  implicit none

! Define output on screen.

  iuOutScr=6
!
! Open input InputFile.
! 
  call getarg(1,InputFile)
!
! Remove the extension .nlam, if present.
!
  if (len_trim(InputFile).gt.5) then
    if (InputFile(len_trim(InputFile)-4:len_trim(InputFile)) &
 &      .eq. '.nlam') InputFile(len_trim(InputFile)-4:)= ' '
  end if
!
  return
 end subroutine shared_header_main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module shared_header

