!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_DAMP. Satish Chimakurthi. 22Aug2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description: Rayleigh damping. This adds a damping matrix. The parameters
!AlfaR, BetaR are assumed to be zero in nlams_input.f90 when they are declared.
!To include this damping, those values can be changed. But, an input card can be created so these
!can be actually read into the code via the input file nlamstest.nlam.
!


module nlams_damp
 contains

 subroutine nlams_dampmat
 use nlams_init, only: K_shell, C_shell, M_shell, K_cg
 use nlams_input, only: AlfaR, BetaR
 implicit none

 C_shell = AlfaR*M_shell + BetaR*K_shell + K_cg

 end subroutine nlams_dampmat

end module nlams_damp
