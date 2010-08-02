!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_solveri_galfa. Satish Chimakurthi. 22Aug2008
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   Time integration of structural dynamic equations using generalized-alfa method
!
!-> Reference:
!
!
!->Subroutines:
!
!     |-nlams_shellsolver_step_galfa()
!
!->Remarks:
!   Generalized-alpha time integration method described in
!   Satish Chimakurthi's dissertation. It was developed based on the Newmark
!   scheme described in the following reference paper:
!
!   Reference paper: A survey of direct time-integration methods in
!   computational structural dynamics - II. Implicit methods.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module nlams_solver_galfa
  real(8),save,allocatable ::DeltaR(:,:)
  real(8),save :: BetaGenAlfa,GammaGenAlfa,alfa_m,alfa_f
  real(8),save :: a0,a1,a2,a3,a4,a5,a6,a7,a8
 
contains
  subroutine nlams_shellsolver_step_galfa()
    use nlams_update
    use nlams_init
    use nlams_input
    use nlams_mass
    use nlams_stiff
    use nlams_damp
    use interface_lapack
    use nlams_interf  
    use mesh_connect_vars, only: NumElDof,Connect,NumDofPerNode,NumNodesPerEl,NumNodes,NumElems,UndefGrid    
    use nlams_loads,       only: bcdof,bcval,NumBC,RHS
!    use eigen
!    use shared_time,       only: FsSolConv
    use nlams_undeformed
    use LociSTREAMdata, only: LociSTREAMRank

    implicit none

! Local variables
    integer::i,j,iErr,kstart,kend,c
    integer,dimension(AnsSize)::AuxLU,AuxLU1
    real(8),dimension(AnsSize,AnsSize)::Keff
    real(8),dimension(AnsSize)::TempDelta
    real(8),dimension(AnsSize,1)::Term1,Term2,Term3,Term4,Term5,Term6,Term7,Term8
    real(8),dimension(:,:),allocatable:: InvM_shell

!Added for convergence check
    integer::CSDConv
    real(8),dimension(AnsSize)::IncrDisp,TotIncrDisp
    real(8)::CSDDelta_pre,norm_curr,norm_pre
    real(8)::CSDDelta,norm1,norm2
    real(8)::norm1_disp,norm1_rot,norm2_disp,norm2_rot
    real(8)::CSDDelta_rot,CSDDelta_disp
    real(8)::NewCSDEps
    real(8)::norm_u,norm_d,CSDDelta_ori

! Energy calculation
!   real(8)::KinetE(1,1),PotenE

! Eigen value analysis
!    real(8),dimension(AnsSize,AnsSize)::AMat

    if(TStep.eq.1) allocate(InvM_shell(AnsSize,AnsSize)); InvM_shell = 0.0

    Keff      = 0.0d0
    TempDelta = 0.0d0
    AuxLU     = 0.0d0
    AuxLU1    = 0.0d0
    
    if(LociSTREAMRank.eq.0) print *, "TimeStep: ", TStep

    if (TStep .gt. 1) then

!---- Preprocessing
! Compute mass, delta_force, stiffness matric based on previous solution
! Compute prescribed motion at previouse time step for recovering data
       call nlams_undeformed_step(Tstep-1)
! M_shell, K_cg, Force_mbd at t-1
       call nlams_massmat

! Recover EleStVec at t-1
          do i=1,NumElems
             do j=1,NumNodesPerEl
                kstart = (Connect(i,j)-1)*NumDofPerNode + 1
                kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
                EleStVec(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=D_curr(kstart:kend,1)
             enddo
          enddo

!       call nlams_update_nodalrotmat

! K_shell, Intforce_g at t-1
       call nlams_stiffmat
! C_shell at t-1
       call nlams_dampmat

! Initialize a new time-step
       call nlams_shellsolver_step_ini()
       call nlams_shellsolver_coeffs

       call nlams_undeformed_step(Tstep)
       call nlams_massmat         
       call nlams_dampmat         
       call nlams_delta_force     
       call nlams_apply_bc        
 
! K_shell_n is the shell stiffness matrix computed at the end of the previous time-step
       Keff = (a0*M_shell) + (a1*C_shell) + (a2*K_shell_n)
       TempDelta = DeltaR(:,1)  

! Calculation of displacement correction
! Solution of the system of equations
       call lapack_lufact(AnsSize,Keff,AuxLU,iErr)
       call lapack_luback(AnsSize,Keff,AuxLU,TempDelta,iErr)

! Displacement increment for 0th iteration of the current time-step 
       Rcorr(:,1) = TempDelta   
 
       Iter = 0
       do while(CSDConv .NE. 1) 
         Iter = Iter + 1
         Dddot_curr(:,1) = (a3*Rcorr(:,1)) - (a4*Ddot_pre(:,1)) - (a5*Dddot_pre(:,1))
         Ddot_curr(:,1)  = (a6*Rcorr(:,1)) - (a7*Ddot_pre(:,1)) - (a8*Dddot_pre(:,1))
         D_curr(:,1) = D_pre(:,1) + Rcorr(:,1)

! Only the displacement correction is stored in the following array
! EleStVecCurr is the vector of global displacements and rotations at all nodes in the structure
! Extracting displacements and rotations of all nodes in an element at the current iteration/time-Iter

         do i=1,NumElems
            do j=1,NumNodesPerEl
               kstart = (Connect(i,j)-1)*NumDofPerNode + 1
               kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
               EleStVecCurr(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=Rcorr(kstart:kend,1)
            enddo
         enddo

! Compute nodal rotation matrices
         call nlams_update_nodalrotmat
 
         FinalAns = D_curr(:,1)

! Total solution at the current time
         DeltaR(:,1) = FinalAns(:)

! Total solution at the current time stored in the following array
! ElemStVec is the vector of global displacements and rotations at all nodes in the structure
! Extracting displacements and rotations of all nodes in an element at the current iteration/time-Iter
          do i=1,NumElems
             do j=1,NumNodesPerEl
                kstart = (Connect(i,j)-1)*NumDofPerNode + 1
                kend   = (Connect(i,j)-1)*NumDofPerNode + NumDofPerNode
                EleStVec(i,(j-1)*NumDofPerNode+1:j*NumDofPerNode)=D_curr(kstart:kend,1)
             enddo
          enddo

          call nlams_stiffmat  !Computes IntForce_g corresponding to the configuration D_curr(:,1)

!Effective stiffness matrix - equation 8/9 of Bonelli's paper
         Keff = (a0*M_shell) + (a1*C_shell) + (a2*K_shell_n)

!Residual force vector  !Equation 10 of Bonelli's paper
!         EqError(:,1) =  ((1.d0 - alfa_f)*Force_mbd(:,1)) + (alfa_f * Force_mbd_pre(:,1)) &
!&                      - ((1.d0 - alfa_m)*matmul(M_shell,Dddot_curr(:,1))) - (alfa_m*matmul(M_shell,Dddot_pre(:,1))) &
!&                      - ((1.d0 - alfa_f)*matmul(C_shell,Ddot_curr(:,1)))  - (alfa_f*matmul(C_shell,Ddot_pre(:,1)))  &
!&                      - ((1.d0 - alfa_f)*IntForce_g(:,1)) - (alfa_f*IntForce_g_n(:,1)) 

          EqError(:,1) = ((1-alfa_f)*AerForce(:,1) + (alfa_f)*AerForce_pre(:,1))+((1-alfa_f)*Force_mbd(:,1)) + &
&                         (alfa_f*Force_mbd_pre(:,1)) - ((1-alfa_m)*matmul(M_shell,Dddot_curr(:,1))) - &
&                         (alfa_m*matmul(M_shell,Dddot_pre(:,1)))-((1-alfa_f)*matmul(C_shell,Ddot_curr(:,1))) - &
&                         (alfa_f*matmul(C_shell,Ddot_pre(:,1))) - ((1-alfa_f)*IntForce_g(:,1)) - (alfa_f*IntForce_g_n(:,1))



!Apply BC to internal force vector
          do i=1,NumBC
            c=bcdof(i)
            EqError(c,1)=bcval(i)
          enddo

          TempDelta = EqError(:,1)   
       
          call lapack_lufact(AnsSize,Keff,AuxLU,iErr)
          if(iErr.eq.1) then
           if (LociSTREAMRank==0)  print  *, 'Error in lapack_lufact'
           stop
          endif

          call lapack_luback(AnsSize,Keff,AuxLU,TempDelta,iErr)
          if(iErr.eq.1) then
           if (LociSTREAMRank==0)  print  *, 'Error in lapack_luback'
           stop
          endif
          
          Rcorr(:,1) = Rcorr(:,1) + TempDelta   !Corrected displacement increment within the time-step

!Added from here for convergence check
          IncrDisp = 0.0
          IncrDisp = TempDelta
          TotIncrDisp = D_pre(:,1) + Rcorr(:,1)

          norm1_disp = 0.0d0
          norm1_rot  = 0.0d0
          norm2_disp = 0.0d0
          norm2_rot  = 0.0d0

          call norm(IncrDisp,AnsSize,norm1_disp,norm1_rot)
          call norm(TotIncrDisp,AnsSize,norm2_disp,norm2_rot)

!          CSDDelta_disp = norm1_disp/ norm2_disp
!          CSDDelta_rot  = norm1_rot / norm2_rot
           CSDDelta_disp = 0.0; CSDDelta_rot = 0.0
           CSDDelta_disp = norm1_disp; CSDDelta_rot = norm1_rot

          call norm_old(IncrDisp,AnsSize,norm_u)
          call norm_old(TotIncrDisp,AnsSize,norm_d)
          
!          CSDDelta     = norm_u/norm_d
          CSDDelta     = 0.0; CSDDelta_ori = 0.0
          CSDDelta_ori = norm_u

        if (LociSTREAMRank==0) then
        write(*,'(i6,1x,i6,4E14.7)')  TStep,Iter,CSDDelta_rot,CSDDelta_disp,CSDDelta,CSDDelta_ori
! File output
          write(2302,'(i6,1xi6,4(E14.7,1x))') TStep,Iter,CSDDelta_rot,CSDDelta_disp,CSDDelta,CSDDelta_ori
        endif
        
!          write(*,'(2i4,3E14.7)') TStep,Iter,CSDDelta_rot,CSDDelta_disp,CSDDelta_pre        
!          if (CSDDelta .LT. CSDMinDelta) then
          if(CSDDelta_ori.le.CSDMinDelta)then
!          if(CSDDelta.lt.CSDMinDelta)then
             if (LociSTREAMRank==0)  print  *,'Solution converged at the current time-step'
             if (LociSTREAMRank==0)  print  *,'Delta at convergence is:',CSDDelta_ori
!              print *,'Delta at convergence is:',CSDDelta
             CSDConv = 1
           
          endif
         if(CSDDelta_ori.eq.1.0)then
           if (LociSTREAMRank==0)  print  *, 'Solution is diverged'
           stop
         endif
!End of addition for convergence check
         if(Iter.ge.5000) then
             if (LociSTREAMRank==0)  print  *, 'It is supposed to be not convergence solution '
             if (LociSTREAMRank==0)  print  *, 'Final Delta X :',CSDDelta_ori
!              print *, 'Final Delta X :',CSDDelta
             stop
          endif
       enddo
       call nlams_deallocate_var
    else
       call nlams_undeformed_step(Tstep)
       call nlams_massmat  
       call nlams_stiffmat
       call nlams_dampmat
       call nlams_delta_force
       call nlams_apply_bc
       call lapack_inv(AnsSize,M_shell,InvM_shell,iErr)
       K_shell_n = K_shell
       if(iErr.eq.1) then
           if (LociSTREAMRank==0)  print  *, 'Error in lapack_inv'
           stop
       endif
!//
        Dddot_curr(:,1) = matmul(InvM_shell,(AerForce(:,1) + Force_mbd(:,1) - matmul(K_shell_n,D_curr(:,1)) - matmul(C_shell,Ddot_curr(:,1))))

!       Dddot_curr(:,1) = matmul(InvM_shell,(Force_mbd(:,1) + RHS(:,1) - matmul(K_shell_n,D_curr(:,1)) - matmul(C_shell,Ddot_curr(:,1))))
! New eigen values and eigenvectors
!       AMat = matmul(InvM_shell,K_shell)
!       call compute_eigen_lapack(AMat,eigenvec,eigenvalue)
!       open(10,file='Eigendata.dat')
!       do i=1,AnsSize
!       write(10,'(i7,1x,f19.9,1x,f19.9)') i,eigenvec(AnsSize-(i-1)),sqrt(eigenvec(AnsSize-(i-1)))
!       enddo
!       close(10)

       call nlams_deallocate_var
!       print *, 'Eignendata has been computed'
!       stop
       deallocate(InvM_shell)
    endif

!Energy calculation
!   KinetE=0.5*matmul(transpose(Ddot_curr),matmul(M_shell,Ddot_curr))
!   PotenE=sum(StrainEn)
!   TotE = KinetE(1,1) + PotenE
!   print *, "Total Energy :",TotE

end subroutine nlams_shellsolver_step_galfa

 subroutine norm(SpanDisp,SpanDispLen,outnorm_disp,outnorm_rot)
  use mesh_connect_vars, only: NumNodes
  implicit none
  integer,intent(in)::SpanDispLen
  real(8),intent(in),dimension(SpanDispLen)::SpanDisp
  integer::i,j
  real(8),intent(out)::outnorm_disp,outnorm_rot

  outnorm_disp = 0.0d0
  outnorm_rot  = 0.0d0

  do i=1,NumNodes
   do j=1,3
   outnorm_disp = outnorm_disp + SpanDisp((i-1)*6+j)*SpanDisp((i-1)*6+j)
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
!
 outnorm = 0.0
!
 do i=1,SpanDispLen
   outnorm = outnorm + SpanDisp(i)**2
 end do
!
 outnorm=sqrt(outnorm)
!
 end subroutine norm_old

  subroutine nlams_shellsolver_coeffs
  use nlams_input, only: DeltaT,SpRadius
  implicit none
  !Generalized-alfa scheme

   alfa_f = ((-SpRadius)/(1+SpRadius))
   alfa_m = ((1-(2*SpRadius))/(1+SpRadius))

!MSC.Marc recommended parameters
!   alfa_f = -0.05
!   alfa_m = 0.0

   BetaGenAlfa  = 0.25d0*(1.d0+alfa_m-alfa_f)**2
   GammaGenAlfa = 0.5d0 + alfa_m - alfa_f

   a0 = (1.d0-alfa_m)/(BetaGenAlfa*DeltaT**2)
   a1 = (((GammaGenAlfa)*(1.d0-alfa_f))/(BetaGenAlfa*DeltaT))
   a2 = (1.d0-alfa_f)
   a3 = 1.d0/(BetaGenAlfa*DeltaT**2)
   a6 = (GammaGenAlfa)/(BetaGenAlfa*DeltaT)
   a4 = 1.d0/(BetaGenAlfa*DeltaT)
   a5 = ((1.d0/(2.d0*BetaGenAlfa)) - 1.d0)
   a7 = ((GammaGenAlfa/BetaGenAlfa)-1.d0)
   a8 = (DeltaT/2.d0)*((GammaGenAlfa/BetaGenAlfa)-2.d0)

  end subroutine nlams_shellsolver_coeffs
 
  subroutine nlams_delta_force
    use mesh_connect_vars
    use nlams_loads, only:RHS
    use nlams_init,  only:IntForce_g_n,Ddot_pre,Dddot_pre,M_shell,C_shell,AnsSize,Force_mbd,K_shell,K_shell_n,Force_mbd_pre,AerForce,AerForce_pre
    use nlams_input, only:num,TStep,DeltaT
!    use interface_lapack
    
    implicit none
    integer :: i,error1
    real(8),dimension(AnsSize,1)::Term1,Term2,Term3,Term4,Term5,Term6,Term7,Term8,Term9,Term10,Term11

    allocate(DeltaR(AnsSize,1)); DeltaR = 0.0

    if (TStep == 1) then
!        DeltaR(:,1) = Force_mbd(:,1) + RHS(:,1) 
        DeltaR(:,1) = Force_mbd(:,1) + AerForce(:,1)
        K_shell_n = K_shell
    else
!Calculation of effective RHS vector
!        Term1(:,1) = ((1-alfa_f)*Force_mbd(:,1)) + ((1-alfa_f)*RHS(:,1))
!        Term2(:,1) = (alfa_f*Force_mbd_pre(:,1)) + (alfa_f*RHS(:,1)) 
!        Term3(:,1) = (((1-alfa_m)/(BetaGenAlfa*DeltaT))*matmul(M_shell,Ddot_pre(:,1)))
!        Term4(:,1) = (((1-alfa_m)/(BetaGenAlfa))*(0.5-BetaGenAlfa))*matmul(M_shell,Dddot_pre(:,1))
!        Term5(:,1) = ((alfa_m)*matmul(M_shell,Dddot_pre(:,1)))
!        Term6(:,1) = ((1-alfa_f)*matmul(C_shell,Ddot_pre(:,1)))
!        Term7(:,1) = ((1-alfa_f)*DeltaT*(1-GammaGenAlfa)*matmul(C_shell,Dddot_pre(:,1)))
!        Term8(:,1) = ((1-alfa_f)*(GammaGenAlfa/BetaGenAlfa)*matmul(C_shell,Ddot_pre(:,1)))
!        Term9(:,1) = ((1-alfa_f)*(GammaGenAlfa*DeltaT/BetaGenAlfa)*(0.5-BetaGenAlfa)*matmul(C_shell,Dddot_pre(:,1)))
!        Term10(:,1) = (1-alfa_f)*IntForce_g_n(:,1)
!        Term11(:,1) = (alfa_f*IntForce_g_n(:,1)) + (alfa_f*matmul(C_shell,Ddot_pre(:,1)))
!        DeltaR(:,1) = Term1(:,1) + Term2(:,1) + Term3(:,1) + Term4(:,1) - Term5(:,1) - Term6(:,1) - Term7(:,1) + Term8(:,1) + Term9(:,1) - Term10(:,1) - Term11(:,1)

        Term1(:,1) = ((1-alfa_f)*Force_mbd(:,1))
        Term2(:,1) = (alfa_f*Force_mbd_pre(:,1))
        Term3(:,1) = (((1-alfa_m)/(BetaGenAlfa*DeltaT))*matmul(M_shell,Ddot_pre(:,1)))
        Term4(:,1) = (((1-alfa_m)/(BetaGenAlfa))*(0.5-BetaGenAlfa))*matmul(M_shell,Dddot_pre(:,1))
        Term5(:,1) = ((alfa_m)*matmul(M_shell,Dddot_pre(:,1)))
        Term6(:,1) = ((1-alfa_f)*matmul(C_shell,Ddot_pre(:,1)))
        Term7(:,1) = ((1-alfa_f)*DeltaT*(1-GammaGenAlfa)*matmul(C_shell,Dddot_pre(:,1)))
        Term8(:,1) = ((1-alfa_f)*(GammaGenAlfa/BetaGenAlfa)*matmul(C_shell,Ddot_pre(:,1)))
        Term9(:,1) = ((1-alfa_f)*(GammaGenAlfa*DeltaT/BetaGenAlfa)*(0.5-BetaGenAlfa)*matmul(C_shell,Dddot_pre(:,1)))
        Term10(:,1) = (1-alfa_f)*IntForce_g_n(:,1)
        Term11(:,1) = (alfa_f*IntForce_g_n(:,1)) + (alfa_f*matmul(C_shell,Ddot_pre(:,1)))
        DeltaR(:,1) = ((1-alfa_f)*AerForce(:,1)) + (alfa_f*AerForce_pre(:,1)) + Term1(:,1) + Term2(:,1) + Term3(:,1) + Term4(:,1) - Term5(:,1) - Term6(:,1) - Term7(:,1) + Term8(:,1) + Term9(:,1) - Term10(:,1) - Term11(:,1)

    endif

  end subroutine nlams_delta_force

  subroutine nlams_shellsolver_step_ini()
    use nlams_init, only: D_pre, D_curr, Ddot_pre, Ddot_curr, Dddot_pre, Dddot_curr, Nodal_CSYS, Nodal_CSYS_pre,PreviousAns,FinalAns,K_shell,K_shell_n,IntForce_g,IntForce_g_n,Force_mbd_pre,Force_mbd,AerForce,AerForce_pre
    implicit none

    D_pre              = D_curr
    Ddot_pre           = Ddot_curr
    Dddot_pre          = Dddot_curr
    K_shell_n          = K_shell
    Nodal_CSYS_pre     = Nodal_CSYS
    PreviousAns        = FinalAns
    IntForce_g_n       = IntForce_g
    Force_mbd_pre(:,1) = Force_mbd(:,1)
!    AerForce_pre(:,1) = AerForce(:,1)

  end subroutine nlams_shellsolver_step_ini
  
! Apply boundary conditions and reduced matrices/vectors
subroutine nlams_apply_bc
    use nlams_init, only: K_shell_n,AnsSize,M_shell,C_shell
    use nlams_loads,only: bcdof,bcval,NumBC
    
    implicit none
    !Local variables
    integer::c,j,i
    
    do i=1,NumBC
       c=bcdof(i)
       do j=1,AnsSize
          K_shell_n(c,j)=0.0
          M_shell(c,j)=0.0
          C_shell(c,j)=0.0
       enddo
       
       K_shell_n(c,c)=1.0
       M_shell(c,c)  =1.0
       C_shell(c,c)  =1.0
       DeltaR(c,1)=bcval(i)
    enddo
      
end subroutine nlams_apply_bc
  
subroutine nlams_deallocate_var
!use nlams_init, only:eigenvec,eigenvlaue
    deallocate(DeltaR)
!    deallocate(eigenvec,eigenvalue)
end subroutine nlams_deallocate_var
  
!subroutine nlams_final_deallocate_var
!    use nlams_loads
!    use mesh_connect_vars
!    deallocate(RHS)
!    deallocate(bcdof)
!    deallocate(bcval)
!    deallocate(Xgl)
!    deallocate(Ygl)
!    deallocate(Zgl)
!    deallocate(UndefGrid)
!    deallocate(Connect)
!end subroutine nlams_final_deallocate_var
  
end module nlams_solver_galfa
