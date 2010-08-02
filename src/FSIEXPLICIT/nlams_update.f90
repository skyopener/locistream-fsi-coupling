module nlams_update
  real(8),save,dimension(3,1)::Xglt_new,Yglt_new,Zglt_new
contains
  
  subroutine nlams_update_local_coord(Xg,Yg,Zg,LocalCSYS)
    use nlams_transforms
    
    implicit none    
    !intent(in) variables
    real(8),dimension(3,1),intent(in)::Xg,Yg,Zg
    
    !Local variables
    real(8)::Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3
    real(8),dimension(3,3)::LocalCSYS
    
    Xg1=Xg(1,1)
    Xg2=Xg(2,1)
    Xg3=Xg(3,1)
    
    Yg1=Yg(1,1)
    Yg2=Yg(2,1)
    Yg3=Yg(3,1)
    
    Zg1=Zg(1,1)
    Zg2=Zg(2,1)
    Zg3=Zg(3,1)
 
    call shell_transf(Xg1,Xg2,Xg3,Yg1,Yg2,Yg3,Zg1,Zg2,Zg3,LocalCSYS)
    
  end subroutine nlams_update_local_coord

! This subroutine is used in nlams_stiff
 subroutine nlams_globalcoord_update(itemp,Xglt1,Yglt1,Zglt1)
    use nlams_init

    implicit none
    !Local variables
    real(8),dimension(3,1)::Ugdisp,Vgdisp,Wgdisp   
    integer :: iElem
    !intent(in) variables
    real(8),dimension(3,1),intent(in)::Xglt1,Yglt1,Zglt1
    integer,intent(in):: itemp

    iElem = 0    
    iElem = itemp
    Ugdisp(1,1) = EleStVec(iElem,1)  !U global displacement of node 1 of element iElem
    Ugdisp(2,1) = EleStVec(iElem,7)  !U global displacement of node 2 of element iElem
    Ugdisp(3,1) = EleStVec(iElem,13) !U global displacement of node 3 of element iElem

    Vgdisp(1,1) = EleStVec(iElem,2)
    Vgdisp(2,1) = EleStVec(iElem,8)
    Vgdisp(3,1) = EleStVec(iElem,14)

    Wgdisp(1,1) = EleStVec(iElem,3)
    Wgdisp(2,1) = EleStVec(iElem,9)
    Wgdisp(3,1) = EleStVec(iElem,15)

    Xglt_new = Xglt1 + Ugdisp
    Yglt_new = Yglt1 + Vgdisp
    Zglt_new = Zglt1 + Wgdisp

  end subroutine nlams_globalcoord_update

! Compute Eq. (1.26) in Satish's SDM paper
! Update nodal rotation matrices
  subroutine nlams_update_nodalrotmat
  use nlams_init, only:Nodal_CSYS,EleStVec,EleStVecCurr,Nodal_CSYS_pre
  use mesh_connect_vars, only:NumElems,NumNodesPerEl
  use nlams_input,only:Iter,TStep
 
  implicit none
  !Local variables
  real(8),dimension(3,3)::Spin_Tens,Ide
  real(8),dimension(3)  ::Pseudo_vec,Rotg
  real(8)::Norm,Denom,Pseudo_vec_norm
  integer::i,j,EleNum,startd

  Ide = 0.0d0
  Ide(1,1) = 1.0d0
  Ide(2,2) = 1.0d0
  Ide(3,3) = 1.0d0
  Denom = 0.d0
  Pseudo_vec_norm = 0.d0
  Spin_Tens = 0.d0

  do EleNum=1,NumElems
    startd = 1
    do i=1,NumNodesPerEl     !i=1,2,3 for a triangle element
!      
      Rotg(1:3) = EleStVecCurr(EleNum,startd+3:startd+5)
! This is incremental rotation of triad S resulted from the last iteration computed in global coordinate system.
!      Norm  =  sqrt(Rotg(1)**2 + Rotg(2)**2 + Rotg(3)**2)
      Pseudo_vec = Rotg(1:3)

! equation (45) in Khosrav et al.'s paper
      Pseudo_vec_norm = sqrt(Pseudo_vec(1)**2 + Pseudo_vec(2)**2 + Pseudo_vec(3)**2)
      Denom = (1.d0 + 0.25d0*Pseudo_vec_norm**2)  
  
! Eq. (1. 29) in Satish's SDM 2009 paper
      Spin_Tens(1,1) =  0.0d0
      Spin_Tens(1,2) = -Pseudo_vec(3)
      Spin_Tens(1,3) =  Pseudo_vec(2)

      Spin_Tens(2,1) =  Pseudo_vec(3)
      Spin_Tens(2,2) =  0.0d0
      Spin_Tens(2,3) = -Pseudo_vec(1)

      Spin_Tens(3,1) = -Pseudo_vec(2)
      Spin_Tens(3,2) =  Pseudo_vec(1)
      Spin_Tens(3,3) =  0.0d0

! Eq. (1.26) in Satish's SDM 2009 paper
      if (TStep .EQ. 1) then
! Matirx Te
        Nodal_CSYS(i,EleNum,:,:) = Ide + (Spin_Tens + 0.5d0*matmul(Spin_Tens,Spin_Tens ))/Denom
      else
! Nodal_CSYS is eq.(2.22) T_bar
        Nodal_CSYS(i,EleNum,:,:) = Ide + (Spin_Tens + 0.5d0*matmul(Spin_Tens,Spin_Tens ))/Denom
! Nodal_CSYS is transformation matrix T_snew
        Nodal_CSYS(i,EleNum,:,:) = matmul(Nodal_CSYS(i,EleNum,:,:),Nodal_CSYS_pre(i,EleNum,:,:))
      endif  
       startd = startd + 6
     enddo
  enddo

 end subroutine nlams_update_nodalrotmat
end module nlams_update



