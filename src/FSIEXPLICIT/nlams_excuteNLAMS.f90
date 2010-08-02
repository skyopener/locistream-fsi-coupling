!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!   Authors:  
!-> Program.- UM/NLAMS Release 1.1
!
!-> Description:  Wing is assumed to be in the X-Y plane.
!------------------------------------------------------------------------------------------------------------------
module nlams_excute
 implicit none

 contains
  subroutine ExcuteNLAMS
 
   use lib_fem       
   use nlams_stiff 
   use nlams_form_shell  !Shell stiffness matrix is formed from DKT and OPT matrices
   use mesh_connect_vars !Reads mesh and connectivity
   use nlams_comp_abbd   !Computes A,B,B,D matrices for composite layups
   use nlams_input       !Contains all subroutines related to NLAMS input
   use nlams_update      !Contains routines to update local element coordinate systems amongst others
   use nlams_init        !Initializes several arrays and allocates memory
   use nlams_solver
   use nlams_solver_galfa
   use nlams_loads
   use nlams_interf
   use nlams_undeformed
   use shared_time, only: shared_time_set
   use nlams_mass
   use nlams_header
   use LociSTREAMdata
   use nlams_outdisp

!  if (LociSTREAMRank==0) print  *, 'before inputLoci'

  call NLAMS_input_LociSTREAM
  if(LociSTREAMRank.eq.0) print  *, 'After nlams_input_LociSTREAM'
  call nlams_header_main
  if(LociSTREAMRank.eq.0) print  *, 'After nlams_header_main'

! Input mesh and connectivity information from LociSTREAM database
  call nlams_mesh_connect
   
  if (LociSTREAMRank==0) print  *, 'After nlams_mesh_connect'
! Compute RHS vector based on bc condition stored LociSTREAM database
  call nlams_nodal_loads
  if (LociSTREAMRank==0) print  *, 'After nlams_nodal_loads'

! Allocation:   
  call nlams_allocate_init(NumNodesPerEl,NumDofPerNode)
  if (LociSTREAMRank==0) print  *, 'After nlams_allocate_inti'

! Compute constitutive matrix:   
  call nlams_abbd              
 
  if (LociSTREAMRank==0) print  *, '**** Start CSD computation ***'

!  if (LociSTREAMRank==0) print  *, 'Maximum Number of Time Step is ', NumTSteps
! Time step loop
!  do TStep = 1,NumTSteps

  if (LociSTREAMRank==0) print  *,'Time step in excuteNLAMS :',TStep

! Set simulation current time
     call shared_time_set (TStep)   

! Set rigid body motions at the current time step
!     call nlams_undeformed_step ()
!     call nlams_undeformed_step

! Integrate forward by a time-step
     if (IntScheme .EQ. 1) then

! Newmark-beta method
       if (LociSTREAMRank==0)  print *,  '*** Newmark-beta method ***'
       call nlams_shellsolver_step ()
     else if (IntScheme .EQ. 2) then

! Generalized-alfa method
       if (LociSTREAMRank==0)  print *,  '*** Generalized-alfa method ***'
       call nlams_shellsolver_step_galfa ()
     
     else
       if (LociSTREAMRank==0) print  *, '*** ERROR in Read IntSchem ***'
       if (LociSTREAMRank==0) print  *, '*** Please check in NLAMS  ***'
       stop
     endif

!  Write output after each time-step
!     if (LociSTREAMRank==0) print *,  '*** a ***', LociSTREAMitfsi
     call nlams_output(FinalAns)
!     if (LociSTREAMRank==0) print *,  '*** b ***', LociSTREAMitfsi
!     print *,  '*** b ***'   
    call NLAMStoLociSTREAM
!     if (LociSTREAMRank==0) print *,  '*** c ***', LociSTREAMitfsi
!     print *,  '*** c ***'
    call nlams_deallocate
!     if (LociSTREAMRank==0) print *,  '*** d ***', LociSTREAMitfsi
!     print *,  '*** d ***'

 end subroutine ExcuteNLAMS
end module nlams_excute       
! end program main

