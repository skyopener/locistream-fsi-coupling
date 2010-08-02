!-> Copyright by The University of Michigan, Aerospace Department. 2002
!
!-> Module.- LIB_KIND
!
!-> Description.-
!
!  This module defines the real(kind) parameters for the compilation.
!  The value returned by the parameter is compiler dependent.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lib_kind
 implicit none
!
 integer,parameter:: real4=selected_real_kind(p=6)
 integer,parameter:: real8=selected_real_kind(p=15)
!
end module lib_kind

