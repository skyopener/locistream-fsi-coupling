!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_solver. Satish Chimakurthi. 22Aug2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   Time integration of structural dynamic equations using Newmark-beta method
!
!-> Reference:
!     Paper: Corotational Nonlinear Analysis of Thin Plates and Shells Using a New Shell Element
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 69, pp. 859-885, 2007
!
!->Subroutines:
!
!     |-nlams_shellsolver_step
!
!->Remarks:
!    Variables are precisely as presented in the reference paper.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_solver
  real(8),save,allocatable ::DeltaR(:,:)
  real(8),save:: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
  
contains
 
 subroutine nlams_shellsolver_step()
    use nlams_update
    use nlams_init
    use nlams_input
    use nlams_mass
    use nlams_stiff
    use nlams_damp
    use interface_lapack
    use nlams_interf  
    use mesh_connect_vars, only: NumElDof,Connect,NumDofPerNode,NumNodesPerEl,NumNodes,NumElems,UndefGrid    
    use nlams_loads,       only: bcdof,bcval,NumBC
    use nlams_undeformed
    use LociSTREAMdata, only: LociSTREAMRank, LociSTREAMitfsi

    implicit none
! Local vaiables
    integer::i,j,iErr,kstart,kend,c,ddj,js,je
    integer,dimension(AnsSize)::AuxLU
    real(8),dimension(AnsSize,AnsSize)::Keff
    real(8),dimension(AnsSize)::TempDelta
    real(8),dimension(AnsSize,1)::EffLoadVec
    real(8),dimension(:,:),allocatable:: InvM_shell


!For convergence check within the Newton-Raphson loop
    real(8),dimension(AnsSize)::IncrDisp,TotIncrDisp
    real(8),dimension(AnsSize)::IniIncrDisp
    real(8)::CSDDelta,norm1_disp,norm1_rot,norm2_disp,norm2_rot
    real(8)::CSDDelta_rot,CSDDelta_disp
    real(8)::NewCSDEps
    real(8)::norm_pre,norm_curr,CSDDelta_pre
    real(8)::norm_u,norm_d,CSDDelta_ori
    integer::CSDConv

!For enery calculation
!    real(8)::KinetE(1,1),PotenE
!    real(8)::KinetE_pre,PotenE_pre,DeltaEnergy

!For LU_Decomp, LU_Solve
!    real(8),dimension(AnsSize)::TempDeltaSol
!    real(8):: max,min,checkvalue
!    real(8),dimension(AnsSize,AnsSize)::Inv_Keff

!    if (LociSTREAMRank==0)  print *,   'TStep :',TStep
    if(TStep.eq.1.and.LociSTREAMitfsi.eq.1) allocate(InvM_shell(AnsSize,AnsSize)); InvM_shell = 0.0

!    print *, "after if(TStep.eq.1)"
 
! Initialization of Local variables
    EffLoadVec  = 0.0d0
    Keff        = 0.0d0
    TempDelta   = 0.0d0
    AuxLU       = 0
 
! TStep is global variable and defined in nlams_input
    if (TStep .GT. 1) then

! ---- Preprocessing 
! Compute mass, delta_force, stiffness matric based on previous solution
! Compute prescribed motion at previouse time step for recovering data
       call nlams_undeformed_step(Tstep-1) 
! M_shell, K_cg, Force_mbd
       call nlams_massmat
! Recover EleStVec:
          do i=1,NumElems
             do j=1,NumNodesPerEl
                kstart = (Connect(i,j)-1)*NumDofPerNode + 1
                kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
                EleStVec(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=D_curr(kstart:kend,1)
             enddo
          enddo
     
!       call nlams_update_nodalrotmat

! K_shell, Intforce_g
       call nlams_stiffmat
!       if (LociSTREAMRank==0)then
!       open(10,file='K_shell.dat')
!       print *, "K_shell:"
!       do i=1,AnsSize
!        print *, K_shell(i,i)
!       enddo
!       close(10)
!       endif 
       
! C_shell
       call nlams_dampmat
! DeltaR
!       call nlams_delta_force

! Initialize a new time-step
       call nlams_shellsolver_step_ini()
       call nlams_shellsolver_coeffs

! Computed reigid body motion
       call nlams_undeformed_step(Tstep)
! M_shell
       call nlams_massmat
       if (LociSTREAMRank==0)  print *,   'after massmat'
!       if (LociSTREAMRank==0)then
!       open(10,file='M_shell.dat')
!       do i=1,AnsSize
!        write(10,*) M_shell(i,i)
!       enddo
!       close(10)
!       endif
     
! C_shell
       call nlams_dampmat
       if (LociSTREAMRank==0)  print *,   'after dampmat'
!       if (LociSTREAMRank==0)then
!       open(10,file='C_shell.dat')
!       do i=1,AnsSize
!        write(10,*) C_shell(i,i)
!       enddo
!       close(10)
 !      endif
       

       
! R vector Step 6 in referece paper
       call nlams_delta_force
       if (LociSTREAMRank==0)  print *,   'after deltaforce'

! Apply boundary condition
       call nlams_apply_bc
       if (LociSTREAMRank==0)  print *,   'after apply_bc'
         
! K_shell_n is the shell stiffness matrix computed at the end of the previous time-step
! Eq. (1.62) in Satish's SDM 2009 paper
       Keff = (a0*M_shell) + (a1*C_shell) + (K_shell_n)
!       if (LociSTREAMRank==0)then
!       open(10,file='K_shell_nbefore.dat')
!         print *, 'K_shell_n:'
!       do i=1,AnsSize
!        !write(10,*) K_shell_n(i,i)
!         print *, K_shell_n(i,i)
!       enddo
!       close(10)
!       endif 
!       if (LociSTREAMRank==0)then
!       open(10,file='K_shell_nbefore.dat')
!         print *, 'K_shell again:'
!       do i=1,AnsSize
!        !write(10,*) K_shell_n(i,i)
!         print *, K_shell(i,i)
!       enddo
!       close(10)
!       endif 
!       Keff = K_shell
!       if (LociSTREAMRank==0)then
!       open(10,file='Keff.dat')
!                print *, 'K_eff:'
!       do i=1,AnsSize
!        print *, Keff(i,i)
!       enddo
!       close(10)
!       endif 

! Eq. (1.63) in Satish's SDM 2009 paper
       TempDelta = DeltaR(:,1)  
       EffLoadVec(:,1) = DeltaR(:,1)

! Calculation of displacement correction, solution of the system of equations	
! Solving Eq.(1.61) in Satish's SDM 2009 paper using LU Decomposition
! LU Decomposition of the global stiffness matrix	
       call lapack_lufact(AnsSize,Keff,AuxLU,iErr)
!        call lu_decomp(AnsSize,Keff,AuxLU,iErr)
       if(iErr.eq.1) then
         if(LociSTREAMRank==0)  print *,  "Error in LAPACK_LUFACT"
         stop
       endif
!        if (LociSTREAMRank==0)  print *,   'after lapack 1'
! LU Back substitution
       call lapack_luback(AnsSize,Keff,AuxLU,TempDelta,iErr)
!        call lu_solve(AnsSize,Keff,AuxLu,TempDelta,TempDeltaSol)      
       if(iErr.eq.1) then
         if (LociSTREAMRank==0)  print *,  "Error in LAPACK_LUBACK"
         stop
       endif
!       if (LociSTREAMRank==0)  print *,   'after lapack 2'
! Displacement increment for 0th iteration of the current time-step - achieved at step 7 of subbaraj's reference paper
! Equation (1.61) in Satish's SDM 2009 paper
       Rcorr(:,1) = TempDelta
       
!       if (LociSTREAMRank==0)then
!         print *, 'Rcorr before NR iteration:'
!         do i=1,AnsSize
!           print *, Rcorr(i,1)
!         enddo
!       endif 
          

! To improve the solution accuracy and to avoid the development of numerical instabilities, it is generally necessary
! to employ iterations within each time step in order to maintain equalibrium.   
       Iter = 0
       do while(CSDConv .NE. 1)
! Step 8(a): i=i+1
          Iter = Iter + 1
! Step 8(b): evaluate the (i-1)th approximation to the acceleration, velocities, and displacements
! Eqs. (1.65) in Satish's SDM 2009 paper
          Dddot_curr(:,1) = (a0*Rcorr(:,1)) - (a2*Ddot_pre(:,1)) - (a3*Dddot_pre(:,1))
          Ddot_curr(:,1)  = (a1*Rcorr(:,1)) - (a4*Ddot_pre(:,1)) - (a5*Dddot_pre(:,1))
          D_curr(:,1) = D_pre(:,1) + Rcorr(:,1)
          
!         if (LociSTREAMRank==0)then
!          print *, 'D_curr, Ddot_curr, Dddot_curr, D_pre:'
!          do i=1,AnsSize
!            print *, D_curr(i,1),Ddot_curr(i,1),Dddot_curr(i,1),D_pre(i,1)
!          enddo
!         endif 
       
! Only the displacement correction is stored in the following array
! EleStVecCurr is the vector of global displacements and rotations at all nodes in the structure
! Extracting displacements and rotations of all nodes in an element at the current iteration/time-Iter
! EleStVecCurr is written in global coordinate system

          do i=1,NumElems
             do j=1,NumNodesPerEl
                kstart = (Connect(i,j)-1)*NumDofPerNode + 1
                kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
                EleStVecCurr(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=Rcorr(kstart:kend,1) 
             enddo
          enddo
                    
! Step b in Satish's SDM 2009 paper: Compute nodal rotation matrices
! Updating Ts, Eq. (44) in reference paper I
          call nlams_update_nodalrotmat
!          if (LociSTREAMRank==0)  print *,   'after update_nodalrotmat'
          
! qt in Eq. (1.69) of Satish's SDM 2009 paper
          FinalAns = D_curr(:,1)
          DeltaR(:,1) = FinalAns(:)
          
!          if (LociSTREAMRank==0)then
!           print *, 'D_curr:'
!           open(10,file='D_curr.dat')
!           do i=1,AnsSize
!             print *, D_curr(i,1)
!              write(10,*) D_curr(i,1)
!           enddo
!           close(10)
!         endif 
          

! Total solution at the current time stored in the following array
! EleStVec is the vector of global displacements and rotations at all nodes in the structure
! Extracting displacements and rotations of all nodes in an element at the current iteration/time-Iter
          do i=1,NumElems
             do j=1,NumNodesPerEl
                kstart = (Connect(i,j)-1)*NumDofPerNode + 1
                kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
                EleStVec(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=D_curr(kstart:kend,1) 
             enddo
          enddo
          
          
          call nlams_stiffmat  !Computes IntForce_g corresponding to the configuration D_curr(:,1)
          
!         if (LociSTREAMRank==0)  print *,   'after stiffmat'
! Residual force vector: Eq. (1.51) in Satish's SDM 2009 paper
! Also Eq.(1.66) in Satish's SDM 2009 paper
!          EqError(:,1) = Force_mbd(:,1)-matmul(M_shell,Dddot_curr(:,1))- matmul(C_shell,Ddot_curr(:,1))-IntForce_g(:,1)
          EqError(:,1) = Force_mbd(:,1) + AerForce(:,1) - matmul(M_shell,Dddot_curr(:,1)) -  matmul(C_shell,Ddot_curr(:,1)) - IntForce_g(:,1)
!           if (LociSTREAMRank==0)  print *,   'after Eq_error 1'
! 
! Apply BC to internal force vector
          do i=1,NumBC
            c=bcdof(i)
            EqError(c,1)=bcval(i)
          enddo
 !                     if (LociSTREAMRank==0)  print *,   'after Eq_error 2'
! TempDelta is RHS vector of Eq. (1.67) in Satish's SDM 2009 paper
          TempDelta = EqError(:,1)
!          if (LociSTREAMRank==0)  print *,   'After Tempdelta'
          Keff = (a0*M_shell) + (a1*C_shell) + (K_shell_n)
!          if (LociSTREAMRank==0)  print *,   'After Keff'

                                

! Corrected displacement increment calculation
!          if (LociSTREAMRank==0)  print *,   'AnsSize before lapack 3: ', AnsSize               
        
          call lapack_lufact(AnsSize,Keff,AuxLU,iErr)
!           if (LociSTREAMRank==0)  print *,   'after lapack 3'
!          call lu_decomp(AnsSize,Keff,AuxLU1,iErr)   
          if(iErr.eq.1) then
            if (LociSTREAMRank==0)  print *,   "Error in LAPACK_LUFACT"
            stop
          endif
          call lapack_luback(AnsSize,Keff,AuxLU,TempDelta,iErr)
!          if (LociSTREAMRank==0)  print *,   'after lapack 4'
!          call lu_solve(AnsSize,Keff,AuxLu1,TempDelta,TempDeltaSol)
          if(iErr.eq.1) then
            if (LociSTREAMRank==0)  print *,   "Error in LAPACK_LUBACK"
            stop
          endif

! Eq.(1.68) in Satish's SDM 2009 paper, evaluate the corrected displacement increment with
          Rcorr(:,1) = Rcorr(:,1) + TempDelta   

          IncrDisp =  0.0
          norm1_disp= 0.0 ;  norm1_rot = 0.0
          norm2_disp= 0.0 ;  norm2_rot = 0.0
          CSDDelta_disp = 0.0; CSDDelta_rot = 0.0; CSDDelta = 0.0; CSDDelta_ori = 0.0
          IncrDisp    = TempDelta
          TotIncrDisp = D_pre(:,1) + Rcorr(:,1)
          
          call norm(IncrDisp,AnsSize,norm1_disp,norm1_rot)
!          if (LociSTREAMRank==0)  print *,   'after norm of rot and disp'
          call norm(TotIncrDisp,AnsSize,norm2_disp,norm2_rot)

!          CSDDelta_disp =(norm1_disp)/(norm2_disp)
!          CSDDelta_rot  =(norm1_rot) /(norm2_rot)
          CSDDelta_disp = norm1_disp; CSDDelta_rot = norm1_rot

          call norm_old(IncrDisp,AnsSize,norm_u)
!          if (LociSTREAMRank==0)  print *,   'after norm of norm_u'
          call norm_old(TotIncrDisp,AnsSize,norm_d)

!          CSDDelta     = norm_u/norm_d
          CSDDelta_ori = norm_u

        if (LociSTREAMRank==0) then
          write(*,'(i6,1x,i6,4E14.6)') TStep,Iter,CSDDelta_disp,CSDDelta_rot,CSDDelta,CSDDelta_ori
          write(2302,'(i6,1x,i6,4(e14.7,1x))')TStep,Iter,CSDDelta_disp,CSDDelta_rot,CSDDelta,CSDDelta_ori
        endif
        
 
          if(CSDDelta_ori.LT.CSDMinDelta)then
!          if(CSDDelta.lt.CSDMinDelta)then
             if (LociSTREAMRank==0)  print *,'Solution converged at the current time-step'
             if (LociSTREAMRank==0)  print *,'Delta at convergence is:',CSDDelta_ori
!              print *,'Delta at convergence is:',CSDDelta
             CSDConv = 1
          endif

!          if(Iter.eq.3) CSDConv = 1

!End of addition for convergence check
          if(Iter.ge.1000) then ! TEMP: Jun 6 2010
             if (LociSTREAMRank==0)  print *, 'It is supposed to be not convergence solution '
             if (LociSTREAMRank==0)  print *, 'Final Delta X :',CSDDelta_ori
!              print *, 'Final Delta X :',CSDDelta
             stop
            ! CSDConv = 1
          endif

          if(CSDDelta_ori.eq.1.0) then
             if (LociSTREAMRank==0)  print *, 'Solution is diverged'
             stop
          endif
       enddo 

       call nlams_deallocate_var
    else if(Tstep.eq.1.and.LociSTREAMitfsi.eq.1)then! Tstep == 1
       call nlams_undeformed_step(Tstep) 
    if (LociSTREAMRank==0)  print *,   'after undeformed_step'
       call nlams_massmat
    if (LociSTREAMRank==0)  print *,   'after massmat'
       call nlams_stiffmat
    if (LociSTREAMRank==0)  print *,   'after stiffmat'
       call nlams_dampmat
    if (LociSTREAMRank==0)  print *,   'after dampmat'
       call nlams_delta_force
    if (LociSTREAMRank==0)  print *,   'after delta force'
       call nlams_apply_bc
    if (LociSTREAMRank==0)  print *,   'after apply_bc'
       call lapack_inv(AnsSize,M_shell,InvM_shell,iErr)
!    if (LociSTREAMRank==0)  print *,   'after lapack_inv'
       if(iErr.eq.1) then
        if (LociSTREAMRank==0)  print *,   'Error in lapack_inv'
        stop
       endif

! Computing accelation of nodes
       Dddot_curr(:,1) = matmul(InvM_shell,(DeltaR(:,1) - matmul(K_shell_n,D_curr(:,1)) - matmul(C_shell,Ddot_curr(:,1))))
       
!       if (LociSTREAMRank==0)then
!        open(10,file='D_Ddot_Dddot.dat')
!        write(10,*) LociSTREAMitfsi
!        do i=1,AnsSize
!       		 write(10,*) D_curr(i,1),Ddot_curr(i,1),Dddot_curr(i,1)
!       	enddo
!       close(10)
!       open(10,file='forcembdaeroforce.dat')
!        do i=1,AnsSize
!       		 write(10,*) Force_mbd(i,1),Aerforce(i,1)
!       	enddo
!       close(10) 
!       endif
       call nlams_deallocate_var
       deallocate(InvM_shell)
    endif
!    if (LociSTREAMRank==0)  print *,   'After solver: ', LociSTREAMitfsi

!Energy calculation
!   KinetE=0.5*matmul(transpose(Ddot_curr),matmul(M_shell,Ddot_curr))
!   PotenE=sum(StrainEn)
!   TotE = KinetE(1,1) + PotenE
!   DeltaEnergy = PotenE - PotenE_pre + KinetE(1,1) - KinetE_pre
!   print *, "Total Energy and DeltaE :",TotE,DeltaEnergy

! Store Kinetic and Pontential Energy
!   KinetE_pre = KinetE(1,1)
!   PotenE_pre = PotenE

end subroutine nlams_shellsolver_step

subroutine norm(SpanDisp,SpanDispLen,outnorm_disp,outnorm_rot)
  use mesh_connect_vars, only: NumNodes
  implicit none

! intent in variables
  integer,intent(in)::SpanDispLen
  real(8),intent(in),dimension(SpanDispLen)::SpanDisp
! intent out variables
  real(8),intent(out)::outnorm_disp,outnorm_rot

! Local variables
  integer:: i,j

  outnorm_disp = 0.0d0
  outnorm_rot  = 0.0d0
  do i=1,NumNodes
   do j=1,3
   outnorm_disp = outnorm_disp + SpanDisp((i-1)*6+j)  *SpanDisp((i-1)*6+j  )
   outnorm_rot  = outnorm_rot  + SpanDisp((i-1)*6+j+3)*SpanDisp((i-1)*6+j+3)
   end do
  enddo

  outnorm_disp=sqrt(outnorm_disp)
  outnorm_rot =sqrt(outnorm_rot)

 end subroutine norm

subroutine norm_old(SpanDisp,SpanDispLen,outnorm)
  implicit none
  integer,intent(in)::SpanDispLen
  real(8),intent(in),dimension(SpanDispLen)::SpanDisp
  integer::i
  real(8),intent(out)::outnorm

  outnorm = 0.0d0
  do i=1,SpanDispLen
    outnorm = outnorm + SpanDisp(i)**2
  end do

  outnorm=sqrt(outnorm)

end subroutine norm_old

! Coefficient for Newmark-Beta method
subroutine nlams_shellsolver_coeffs
  use nlams_input, only: BetaNewmark, GammaNewmark, DeltaT
  implicit none

   a0 = 0.0d0; a1 = 0.0d0; a2 = 0.0d0; a3= 0.0d0
   a4 = 0.0d0; a5 = 0.0d0; a6 = 0.0d0; a7= 0.0d0
   a8 = 0.0d0; a9 = 0.0d0; a10= 0.0d0

   a0  = 1.d0/(BetaNewmark*DeltaT*DeltaT)
   a1  = GammaNewmark/(BetaNewmark*DeltaT)
   a2  = 1.d0/(BetaNewmark*DeltaT)
   a3  = ((1.d0/(2.d0*BetaNewmark)) - 1.d0)
   a4  = ((GammaNewmark/BetaNewmark)-1.d0)
   a5  = (DeltaT/2.d0)*((GammaNewmark/BetaNewmark)-2.d0)
   a6  =  a0
   a7  = -a2
   a8  = -a3
   a9  = (DeltaT*(1.d0-GammaNewmark))
   a10 = (DeltaT*GammaNewmark)

end subroutine nlams_shellsolver_coeffs

! Store the data computed at current timestep 
subroutine nlams_shellsolver_step_ini()
   use LociSTREAMdata, only: LociSTREAMRank
   use nlams_init, only: AnsSize, D_pre, D_curr, Ddot_pre, Ddot_curr, Dddot_pre, Dddot_curr, Nodal_CSYS, Nodal_CSYS_pre,PreviousAns,FinalAns,K_shell,K_shell_n,IntForce_g,IntForce_g_n
   implicit none
   integer :: i
 
    D_pre          = D_curr    ! Displacement state vector
    Ddot_pre       = Ddot_curr ! Velocity state vector
    Dddot_pre      = Dddot_curr! Acceleration state vector
    K_shell_n      = K_shell   ! Shell stiffness matrix
    Nodal_CSYS_pre = Nodal_CSYS ! Nodal coordinate system
    PreviousAns    = FinalAns ! Total state vector of the current time step
    IntForce_g_n   = IntForce_g

end subroutine nlams_shellsolver_step_ini
  
! Compute load inbalance 
  subroutine nlams_delta_force 
    use mesh_connect_vars
!    use nlams_loads, only:RHS
    use nlams_init,  only:IntForce_g_n,Ddot_pre,Dddot_pre,M_shell,C_shell,AnsSize,Force_mbd,K_shell,K_shell_n,AerForce, AnsSize
    use nlams_input, only:TStep,DeltaT
    use LociSTREAMdata, only: LociSTREAMRank

    implicit none    
    integer :: i,error1
    real(8),dimension(AnsSize,1)::Term1,Term2
    
    allocate(DeltaR(AnsSize,1)); DeltaR = 0.0d0
   
!    print *, 'TStep: ',TStep 
    if (TStep == 1) then
!        DeltaR(:,1) = Force_mbd(:,1) 
        DeltaR(:,1) = Force_mbd(:,1) + AerForce(:,1)
        K_shell_n   = K_shell
    else
!    if(LociSTREAMRank==0) print *, 'inside deltaR'
        Term1(:,1) = a2*Ddot_pre(:,1) + a3*Dddot_pre(:,1)
        Term2(:,1) = a4*Ddot_pre(:,1) + a5*Dddot_pre(:,1)
!        DeltaR(:,1)= Force_mbd(:,1) + matmul(M_shell,Term1(:,1)) + matmul(C_shell,Term2(:,1)) - IntForce_g_n(:,1)
        DeltaR(:,1)= Force_mbd(:,1) + AerForce(:,1) + matmul(M_shell,Term1(:,1)) + matmul(C_shell,Term2(:,1)) - IntForce_g_n(:,1)

    endif
 
  end subroutine nlams_delta_force
  
! Apply boundary conditions and reduced matrices/vectors
  subroutine nlams_apply_bc
    use nlams_init, only: K_shell_n,AnsSize,M_shell,C_shell
    use nlams_loads,only: bcdof,bcval,NumBC
!    use nlams_input,only: TStep
    implicit none    

!Local variables
    integer::c,i,j
    c = 0
   
    do i=1,NumBC
       c=bcdof(i)
       do j=1,AnsSize
          K_shell_n(c,j)= 0.0d0
          M_shell(c,j)  = 0.0d0
          C_shell(c,j)  = 0.0d0
       enddo
 
       K_shell_n(c,c)=1.0d0
       M_shell(c,c)  =1.0d0
       C_shell(c,c)  =1.0d0
       DeltaR(c,1)=bcval(i)
    enddo
     
 end subroutine nlams_apply_bc

! Deallocation of DeltaR  
  subroutine nlams_deallocate_var
    implicit none
    deallocate(DeltaR);
  end subroutine nlams_deallocate_var
  
! Deallocation of RHS,bcdof,bcval,Xgl,Ygl,Zgl,UndefGrid,Connect
! subroutine nlams_final_deallocate_var
!    use nlams_loads
!    use mesh_connect_vars
!    implicit none
!    deallocate(RHS)
!    deallocate(bcdof)
!    deallocate(bcval)
!    deallocate(Xgl)
!    deallocate(Ygl)
!    deallocate(Zgl)
!    deallocate(UndefGrid)
!    deallocate(Connect)
! end subroutine nlams_final_deallocate_var
 
end module nlams_solver
