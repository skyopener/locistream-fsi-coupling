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

    implicit none
! Local vaiables
    integer::i,j,iErr,kstart,kend,c,ddj,js,je
    integer,dimension(AnsSize)::AuxLU
    real(8),dimension(AnsSize,AnsSize)::Keff
    real(8),dimension(AnsSize)::TempDelta
!    real(8),dimension(AnsSize,1)::EffLoadVec
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

!    if(TStep.eq.1) allocate(InvM_shell(AnsSize,AnsSize)); InvM_shell = 0.0

!    print *, "after if(TStep.eq.1)"
 
! Initialization of Local variables
!    EffLoadVec  = 0.0d0
    Keff        = 0.0d0
    TempDelta   = 0.0d0
    AuxLU       = 0.0d0
 
! TStep is global variable and defined in nlams_input
!    if (TStep .GT. 1) then
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Compute mass, delta_force, stiffness matric based on previous solution
! Compute prescribed motion at previouse time step for recovering data
 if(TStep.gt.1)then
       call nlams_undeformed_step(Tstep-1)  ! Rigid body motion at (it-1)
       call nlams_massmat                   ! M_shell, K_cg, Force_mbd at (it-1)

       do i=1,NumElems
        do j=1,NumNodesPerEl
          kstart = (Connect(i,j)-1)*NumDofPerNode + 1
          kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
          EleStVec(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=D_curr(kstart:kend,1)
        enddo
       enddo

       call nlams_stiffmat ! K_shell, Intforce_g
       call nlams_dampmat  ! C_shell
! DeltaR
!     call nlams_delta_force
 endif
       call nlams_shellsolver_step_ini() ! Updating D_pre, Ddot_pre, Dddot_pre, k_shell_n, Nodal_CSYS_pre, PreviousAns, IntForce_g
       call nlams_shellsolver_coeffs     ! 

       call nlams_undeformed_step(Tstep) ! Rigid body motion at it
       call nlams_massmat                ! M_shell, k_cg, and Force_mbd at it
       call nlams_dampmat                ! C_shell at it
       call nlams_delta_force            ! DetalR at it 
       call nlams_apply_bc               ! Apply boundary condition
       
       Keff = (a0*M_shell) + (a1*C_shell) + (K_shell_n) ! Effective stiffness matrix at it
       TempDelta = DeltaR(:,1)                          ! Effective load vector at it
 
!      EffLoadVec(:,1) = DeltaR(:,1) ! Effective

! Solving Equation (2.56) in Satish's thesis
       call lapack_lufact(AnsSize,Keff,AuxLU,iErr) ! LU Decomposition
!        call lu_decomp(AnsSize,Keff,AuxLU,iErr)
       if(iErr.eq.1) then
         write(*,*) "Error in LAPACK_LUFACT"
         stop
       endif
       call lapack_luback(AnsSize,Keff,AuxLU,TempDelta,iErr) ! LU Backword 
!        call lu_solve(AnsSize,Keff,AuxLu,TempDelta,TempDeltaSol)      
       if(iErr.eq.1) then
         write(*,*) "Error in LAPACK_LUBACK"
         stop
       endif

       Rcorr(:,1) = TempDelta ! Rcorr is genelized displacement vecor at it & nr=0

     
! Until here time level loop 
!----------------------------------------------------------------------------------------
! Since here each iteralation level loop
!       Iter = 0
!       do while(CSDConv .NE. 1)
!       Iter = Iter + 1

! Evaluate the (i-1)th approximation to the acceleration, velocities, and displacements
        Dddot_curr(:,1) = (a0*Rcorr(:,1)) - (a2*Ddot_pre(:,1)) - (a3*Dddot_pre(:,1))
        Ddot_curr(:,1)  = (a1*Rcorr(:,1)) - (a4*Ddot_pre(:,1)) - (a5*Dddot_pre(:,1))
        D_curr(:,1)     = D_pre(:,1) + Rcorr(:,1)
       
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
                    
      call nlams_update_nodalrotmat ! Updating Ts, Eq. (44) in reference paper I
          
      FinalAns = D_curr(:,1)  
      DeltaR(:,1) = FinalAns(:)

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
         
! Residual force vector: Eq. (2.61) in Satish's thesis
! 1) No aerodynamic force
!     EqError(:,1) = Force_mbd(:,1)-matmul(M_shell,Dddot_curr(:,1))- matmul(C_shell,Ddot_curr(:,1))-IntForce_g(:,1)
! 2) With aerodynamic force 
! In Loci-STREAM and NLAMS coupling, aeroforce(:,1) is changing iteration level, however, Force_mbd is not changing in iteration level
      EqError(:,1) = Force_mbd(:,1) + AerForce(:,1) - matmul(M_shell,Dddot_curr(:,1)) -  matmul(C_shell,Ddot_curr(:,1)) - IntForce_g(:,1)
 
! Apply BC to internal force vector
      do i=1,NumBC
        c=bcdof(i)
        EqError(c,1)=bcval(i)
      enddo
           
      TempDelta = EqError(:,1)
      Keff = (a0*M_shell) + (a1*C_shell) + (K_shell_n)

! Corrected displacement increment calculation, solving Equation (2.62) in Satish's thesis
          call lapack_lufact(AnsSize,Keff,AuxLU,iErr)
!          call lu_decomp(AnsSize,Keff,AuxLU1,iErr)   
          if(iErr.eq.1) then
            write(*,*) "Error in LAPACK_LUFACT"
            stop
          endif

          call lapack_luback(AnsSize,Keff,AuxLU,TempDelta,iErr)
!          call lu_solve(AnsSize,Keff,AuxLu1,TempDelta,TempDeltaSol)
          if(iErr.eq.1) then
            write(*,*) "Error in LAPACK_LUBACK"
            stop
          endif

! Eq.(1.68) in Satish's SDM 2009 paper, evaluate the corrected displacement increment with
          Rcorr(:,1) = Rcorr(:,1) + TempDelta   

          IncrDisp =  0.0
          norm1_disp= 0.0 ;  norm1_rot = 0.0
          norm2_disp= 0.0 ;  norm2_rot = 0.0

          IncrDisp    = TempDelta
          TotIncrDisp = D_pre(:,1) + Rcorr(:,1)
          
          call norm(IncrDisp,AnsSize,norm1_disp,norm1_rot)
          call norm(TotIncrDisp,AnsSize,norm2_disp,norm2_rot)

          CSDDelta_disp =(norm1_disp)/(norm2_disp)
          CSDDelta_rot  =(norm1_rot) /(norm2_rot)

          call norm_old(IncrDisp,AnsSize,norm_u)
          call norm_old(TotIncrDisp,AnsSize,norm_d)

          CSDDelta     = norm_u/norm_d
          CSDDelta_ori = norm_u

          write(*,'(i6,1x,i6,4E14.6)') TStep,Iter,CSDDelta_disp,CSDDelta_rot,CSDDelta,CSDDelta_ori
!          write(2302,'(i6,1x,i6,4(e14.7,1x))')TStep,Iter,CSDDelta_disp,CSDDelta_rot,CSDDelta,CSDDelta_ori
!          if(CSDDelta_ori.LT.CSDMinDelta)then
!          if(CSDDelta.lt.CSDMinDelta)then
!             print *,'Solution converged at the current time-step'
!             print *,'Delta at convergence is:',CSDDelta_ori
!              print *,'Delta at convergence is:',CSDDelta
!             CSDConv = 1
!          endif

!End of addition for convergence check
!          if(Iter.ge.5000) then
!             print *, 'It is supposed to be not convergence solution '
!             print *, 'Final Delta X :',CSDDelta_ori
!              print *, 'Final Delta X :',CSDDelta
!             stop
!          endif

! Error: 
          if(CSDDelta.ge.1.0) then
             print *, 'Solution is diverged'
             stop
          endif
!       enddo 

       call nlams_deallocate_var

!    else! Tstep == 1
!       call nlams_undeformed_step(Tstep) 
!    write(*,*) 'after undeformed_step'
!       call nlams_massmat
!    write(*,*) 'after massmat'
!       call nlams_stiffmat
!    write(*,*) 'after stiffmat'
!       call nlams_dampmat
!    write(*,*) 'after dampmat'
!       call nlams_delta_force
!    write(*,*) 'after delta force'
!       call nlams_apply_bc
!    write(*,*) 'after apply_bc'
!       call lapack_inv(AnsSize,M_shell,InvM_shell,iErr)
!       if(iErr.eq.1) then
!        write(*,*) 'Error in lapack_inv'
!        stop
!       endif
! Computing accelation of nodes
!       Dddot_curr(:,1) = matmul(InvM_shell,(DeltaR(:,1) - matmul(K_shell_n,D_curr(:,1)) - matmul(C_shell,Ddot_curr(:,1))))
!
!       open(10,file='M_Shell.dat')
!       do i=1,AnsSize
!        write(10,*) M_Shell(i,1)
!       enddo
!       close(10)
!       open(10,file='Dddot_curr.dat')
!       do i=1,AnsSize
!        write(10,*) Dddot_curr(i,1),Ddot_curr(i,1),D_curr(i,1)
!       enddo
!       close(10)
!
!
!       call nlams_deallocate_var
!       deallocate(InvM_shell)
!    endif

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
   use nlams_init, only: D_pre, D_curr, Ddot_pre, Ddot_curr, Dddot_pre, Dddot_curr, Nodal_CSYS, Nodal_CSYS_pre,PreviousAns,FinalAns,K_shell,K_shell_n,IntForce_g,IntForce_g_n
   implicit none
 
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

    implicit none    
    integer :: i,error1
    real(8),dimension(AnsSize,1)::Term1,Term2
    
    allocate(DeltaR(AnsSize,1)); DeltaR = 0.0d0

    Term1 = 0.d0; Term2 = 0.d0
 
    Term1(:,1) = a2*Ddot_pre(:,1) + a3*Dddot_pre(:,1)
    Term2(:,1) = a4*Ddot_pre(:,1) + a5*Dddot_pre(:,1)   
    DeltaR(:,1)= Force_mbd(:,1) + AerForce(:,1) + matmul(M_shell,Term1(:,1)) + matmul(C_shell,Term2(:,1)) - IntForce_g_n(:,1)
   
!    print *, 'TStep: ',TStep 
!    if (TStep == 1) then
!        DeltaR(:,1) = Force_mbd(:,1) 
!        DeltaR(:,1) = Force_mbd(:,1) + AerForce(:,1)
!        K_shell_n   = K_shell
!    else
!        Term1(:,1) = a2*Ddot_pre(:,1) + a3*Dddot_pre(:,1)
!        Term2(:,1) = a4*Ddot_pre(:,1) + a5*Dddot_pre(:,1)
!        DeltaR(:,1)= Force_mbd(:,1) + matmul(M_shell,Term1(:,1)) + matmul(C_shell,Term2(:,1)) - IntForce_g_n(:,1)
!        DeltaR(:,1)= Force_mbd(:,1) + AerForce(:,1) + matmul(M_shell,Term1(:,1)) + matmul(C_shell,Term2(:,1)) - IntForce_g_n(:,1)
!    endif
 
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
