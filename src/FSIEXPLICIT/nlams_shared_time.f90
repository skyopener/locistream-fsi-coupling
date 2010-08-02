!-> Copyright by The University of Michigan, Aerospace Department. 2004
!
!-> Module.- SHARED_TIME. Rafa Palacios. 28Jul2004
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!    Initialize different modules in the program.
!
!-> Contains.-
!
!    shared_time_input:  Read input time from file.
!    shared_time_set:    Set current time step.
!    shared_time_feval:  Evaluate time function.
!    shared_time_fdiff:  Evaluate differential of time function.
!    shared_time_finteg: Evaluate integral of time function.

! Satish:
!    shared_time_cordrot_flapfourier:  Fourier series parameters for prescribed flap motion
!    shared_time_cordrot_lagfourier:   Fourier series parameters for prescribed lag motion
!    shared_time_cordrot_pitchfourier: Fourier series parameters for pitching motion

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module shared_time
  use mod_shared
  implicit none

 real(8),save:: TheTime                         ! Current time of the simulation.

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!-> Subroutine SHARED_TIME_SET
!
!-> Description.- 
! 
!     Set the current time step and total time for the current step. It also
!     stores time at previous steps and sets coefficients for finite-difference
!     approximations to time derivatives.
!
!-> Remarks.-
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine shared_time_set (iStep)
  use nlams_input

! Definition of variables.

  integer,intent(in):: iStep     ! Current step.

  TheTime =  iStep   *DeltaT
  
  return
 end subroutine shared_time_set
end module shared_time
