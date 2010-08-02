!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_undeformed. Satish Chimakurthi. 2Jan2009
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   The Cartesian rotation vector is defined by as the vector which has the direction
!   of the rotation axis and a length equal to the amplitude of the rotation.  
!   The rotational operator can be expressed directly either in trigonometric form starting
!   from the general expression 
!   T = I + sin(norm of fai)/norm of fai*skew(fai) + 1-cos(norm of fai)/(norm of fai)^2*skew(fai)*skew(fai)
!   Transpose T will transform the components of a vector in global to inertial frame
!-> Reference:
!
!->Subroutines:
!
!     |-nlams_undeformed_step
!
!
!->Remarks:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nlams_undeformed
real(8),save,dimension(3,3)::T,SkewOmega,SkewOmegadot,tempval
real(8),save,dimension(18,18)::Tf
real(8),save,dimension(3,1)::XRdot,MemOri

contains
  
 subroutine nlams_undeformed_step(it)
 use mesh_connect_vars, only: NumElDof
 use nlams_input, only: DeltaT,TStep,PitchAmp,FlapAmp,LagAmp,FlapOmega,RotType,PlungeOmega,XfiAmp,YfiAmp,ZfiAmp,PlugType
! use mod_shared,        only: Pi
 use lib_rodrigues
 use LociSTREAMdata, only: LociSTREAMRank, LociSTREAMitfsi
! use shared_time,       only: time

 !Local variables
 real(8),dimension(3)::RotVec,RotVecdot,RotVecddot
 real(8),dimension(3)::RotVec_pre,RotVecdot_pre,RotVecddot_pre
 real(8),dimension(3)::RotVec_cur,RotVecdot_cur,RotVecddot_cur

 integer, intent(in) :: it

 real(8) :: time, Omega

 integer ::k, i
! real(8)::sin

!For plunge motion
 real(8)::Xfi,Yfi,Zfi,Xfidot,Yfidot,Zfidot
 
!Delay coefficient
 real(8):: value
 real(8):: Pi

 Pi = acos(-1.0)
 time = DeltaT * real(it)
 Xfi = 0.d0; Yfi = 0.d0; Zfi = 0.d0
 Xfidot = 0.d0 ; Yfidot = 0.d0 ; Zfidot = 0.d0

! Evaluate undeformed body kinematics (rotation vector at the current time-step and transformation matrix)  
! 1-axis is along the length of the wing
! 2-axis is along the width 
! 3-axis is normal to the wing surface

! Coefficients
! time is current computational time
  Omega = FlapOmega*2.*Pi
  value = 4000
  
!Omega is the frequency of excitation, it is the same for all three DOF i.e. Pitch, flap, and lag
  RotVec    = 0.0d0
  RotVecdot = 0.0d0
  RotVecddot= 0.0d0
 
!  print *, sin(Omega*time),cos(Omega*time),Omega**2

! NEW
!  if(time.ge.(2.0/FlapOmega))then
!    Omega = 0.d0
!endif 
!  print *, 'FlapAmp in undeformed: ', FlapAmp
  if(LociSTREAMRank==0) print *, 'PlugType, RotType: ',PlugType, RotType
  if(LociSTREAMRank==0) print *, 'PitchAmp, FlapAmp, LagAmp: ',PitchAmp, FlapAmp, LagAmp
  if(LociSTREAMRank==0) print *, 'XfiAmp,YfiAmp,ZfiAmp: ',XfiAmp,YfiAmp,ZfiAmp

!  if(LociSTREAMRank==0) read(*,*)

  if (PlugType.eq.0.and.RotType.ne.0)then
  if (RotType == 1) then  !sine excitation
!  if(LociSTREAMRank==0) print *, 'RotType= 1'
  RotVec(1)    = (PitchAmp*Pi/180.0)*sin(Omega*time) !Current pitch angle (in radians)
  RotVecdot(1) = (PitchAmp*Pi/180.0)*Omega*cos(Omega*time)
  RotVecddot(1)=-(PitchAmp*Pi/180.0)*(Omega**2)*sin(Omega*time)
 
  RotVec(2)    = (FlapAmp*Pi/180.0)*sin(Omega*time)  !Current flap angle (in radians)
  RotVecdot(2) = (FlapAmp*Pi/180.0)*Omega*cos(Omega*time)
  RotVecddot(2)=-(FlapAmp*Pi/180.0)*(Omega**2)*sin(Omega*time)
 
  RotVec(3)    = (LagAmp*Pi/180.0)*sin(Omega*time)   !Current lag angle (in radians)
  RotVecdot(3) = (LagAmp*Pi/180.0)*Omega*cos(Omega*time)
  RotVecddot(3)=-(LagAmp*Pi/180.0)*(Omega**2)*sin(Omega*time)

 else if(RotType==3)then  !1-cos excitation
!  if(LociSTREAMRank==0) print *, 'RotType= 3'
  RotVec(1)    =(PitchAmp*Pi/180.0)*(1.0 - cos(Omega*time))                    !Current flap angle (in radians)
  RotVecdot(1) =(PitchAmp*Pi/180.0)*Omega*sin(Omega*time)
  RotVecddot(1)=(PitchAmp*Pi/180.0)*(Omega**2)*cos(Omega*time)

  RotVec(2)    =(FlapAmp*Pi/180.0)*(1.0 - cos(Omega*time))                    !Current flap angle (in radians)
  RotVecdot(2) =(FlapAmp*Pi/180.0)*Omega*sin(Omega*time)
  RotVecddot(2)=(FlapAmp*Pi/180.0)*(Omega**2)*cos(Omega*time)

  RotVec(3)    =(LagAmp*Pi/180.0)*(1.0 - cos(Omega*time))                    !Current flap angle (in radians)
  RotVecdot(3) =(LagAmp*Pi/180.0)*Omega*sin(Omega*time)
  RotVecddot(3)=(LagAmp*Pi/180.0)*(Omega**2)*cos(Omega*time)

 else if(RotType==2)then
!  if(LociSTREAMRank==0) print *, 'RotType= 2'
! cosine excitaion
  RotVec(1)    = (PitchAmp*Pi/180.0)*(cos(Omega*time))                    !Current flap angle (in radians)
  RotVecdot(1) =-(PitchAmp*Pi/180.0)*Omega*sin(Omega*time)
  RotVecddot(1)=-(PitchAmp*Pi/180.0)*(Omega**2)*cos(Omega*time)

  RotVec(2)    = (FlapAmp*Pi/180.0)*cos(Omega*time)                    !Current flap angle (in radians)
  RotVecdot(2) =-(FlapAmp*Pi/180.0)*Omega*sin(Omega*time)
  RotVecddot(2)=-(FlapAmp*Pi/180.0)*(Omega**2)*cos(Omega*time)

  RotVec(3)    = (LagAmp*Pi/180.0)*cos(Omega*time)                   !Current flap angle (in radians)
  RotVecdot(3) =-(LagAmp*Pi/180.0)*Omega*sin(Omega*time)
  RotVecddot(3)=-(LagAmp*Pi/180.0)*(Omega**2)*cos(Omega*time)

 else if(RotType==4)then
! (1-exp(-value*t^2))*sin excitation
! 
  RotVec(1) = (PitchAmp*Pi/180)*(1-exp(-value*time*time))*sin(Omega*time)
  RotVecdot(1) = ((PitchAmp*Pi/180)*(2*time*value)*(exp(-value*time*time))*sin(Omega*time)) +&
&  ((PitchAmp*Pi/180)*(Omega)*cos(Omega*time)*(1-exp(-value*time*time)))
  RotVecddot(1) = ((PitchAmp*Pi/180)*(2*value)*(exp(-value*time*time))*sin(Omega*time)) + &
& ((PitchAmp*Pi/180)*(2*time*value)*(-2*time*value)*(exp(-value*time*time))*sin(Omega*time)) &
& + ((PitchAmp*Pi/180)*(2*time*value)*(exp(-value*time*time))*(Omega)*cos(Omega*time)) &
& - ((PitchAmp*Pi/180)*(Omega)*(Omega)*sin(Omega*time)*(1-exp(-value*time*time))) &
& + ((PitchAmp*Pi/180)*(Omega)*cos(Omega*time)*(-2*value*time)*(1-exp(-value*time*time)))

  RotVec(2) = (FlapAmp*Pi/180)*(1-exp(-value*time*time))*sin(Omega*time)
  RotVecdot(2) = ((FlapAmp*Pi/180)*(2*time*value)*(exp(-value*time*time))*sin(Omega*time)) &
& +  ((FlapAmp*Pi/180)*(Omega)*cos(Omega*time)*(1-exp(-value*time*time)))
  RotVecddot(2) = ((FlapAmp*Pi/180)*(2*value)*(exp(-value*time*time))*sin(Omega*time)) + &
& ((FlapAmp*Pi/180)*(2*time*value)*(-2*time*value)*(exp(-value*time*time))*sin(Omega*time)) &
& + ((FlapAmp*Pi/180)*(2*time*value)*(exp(-value*time*time))*(Omega)*cos(Omega*time)) &
& - ((FlapAmp*Pi/180)*(Omega)*(Omega)*sin(Omega*time)*(1-exp(-value*time*time))) &
& + ((FlapAmp*Pi/180)*(Omega)*cos(Omega*time)*(-2*value*time)*(1-exp(-value*time*time)))

  RotVec(3) = (LagAmp*Pi/180)*(1-exp(-value*time*time))*sin(Omega*time)
  RotVecdot(3) = ((LagAmp*Pi/180)*(2*time*value)*(exp(-value*time*time))*sin(Omega*time)) &
& + ((LagAmp*Pi/180)*(Omega)*cos(Omega*time)*(1-exp(-value*time*time)))
  RotVecddot(3) = ((LagAmp*Pi/180)*(2*value)*(exp(-value*time*time))*sin(Omega*time)) + &
& ((LagAmp*Pi/180)*(2*time*value)*(-2*time*value)*(exp(-value*time*time))*sin(Omega*time)) &
& + ((LagAmp*Pi/180)*(2*time*value)*(exp(-value*time*time))*(Omega)*cos(Omega*time)) &
& - ((LagAmp*Pi/180)*(Omega)*(Omega)*sin(Omega*time)*(1-exp(-value*time*time))) &
& + ((LagAmp*Pi/180)*(Omega)*cos(Omega*time)*(-2*value*time)*(1-exp(-value*time*time)))
 endif
endif

! 3x3 dynamic transformation matrix 
! using function rodr_rvec2mat
   T = rodr_rvec2mat(RotVec)

   if (LociSTREAMRank==0)  print *, 'Flap/Plung Omega, time : ', Omega, time   

!   do i=1,3
!     write(*,*) T(i,1:3)
!   enddo
!   read(*,*)

!Transformation matrix from global (flapping) to inertial frame
    k = 1
    do i=1,6
       Tf(k:k+2,k:k+2) = T(1:3,1:3)
       k = k + 3
    enddo

!   do i=1,18
!     write(*,*) Tf(i,1:18)
!   enddo
!   read(*,*)

  SkewOmega   =vect_skew(-RotVecdot)
  SkewOmegadot=vect_skew(-RotVecddot)

!   do i=1,3
!     write(*,*) SkewOmega(i,1:3)
!   enddo
!   read(*,*)
!      do i=1,3
!     write(*,*) SkewOmegadot(i,1:3)
!   enddo
!   read(*,*)


!Evaluate rigid body translational velocities and accelerations based on prescribed position of the origin of
!flapping wing coordinate system w.r.t. the inertial frame

 if(RotType.eq.0.and.PlugType.ne.0)then
 if(PlugType.eq.3)then
 Xfi = (XfiAmp)*(1-cos(2*Pi*PlungeOmega*time))    
 Yfi = (YfiAmp)*(1-cos(2*Pi*PlungeOmega*time))
 Zfi = (ZfiAmp)*(1-cos(2*Pi*PlungeOmega*time))

 MemOri(1,1) = Xfi
 MemOri(2,1) = Yfi
 MemOri(3,1) = Zfi

 Xfidot = (XfiAmp)*2*Pi*PlungeOmega*sin(2*Pi*PlungeOmega*time)
 Yfidot = (YfiAmp)*2*Pi*PlungeOmega*sin(2*Pi*PlungeOmega*time)
 Zfidot = (ZfiAmp)*2*Pi*PlungeOmega*sin(2*Pi*PlungeOmega*time)

 XRdot(1,1)=(XfiAmp)*4*(Pi**2)*(PlungeOmega**2)*cos(2*Pi*PlungeOmega*time)
 XRdot(2,1)=(YfiAmp)*4*(Pi**2)*(PlungeOmega**2)*cos(2*Pi*PlungeOmega*time)
 XRdot(3,1)=(ZfiAmp)*4*(Pi**2)*(PlungeOmega**2)*cos(2*Pi*PlungeOmega*time)
 elseif(PlugType.eq.1)then
 Xfi = (XfiAmp)*sin(2*Pi*PlungeOmega*time)
 Yfi = (YfiAmp)*sin(2*Pi*PlungeOmega*time)
 Zfi = (ZfiAmp)*sin(2*Pi*PlungeOmega*time)

 MemOri(1,1) = Xfi
 MemOri(2,1) = Yfi
 MemOri(3,1) = Zfi

 Xfidot = (XfiAmp)*2.*Pi*PlungeOmega*cos(2*Pi*PlungeOmega*time)
 Yfidot = (YfiAmp)*2.*Pi*PlungeOmega*cos(2*Pi*PlungeOmega*time)
 Zfidot = (ZfiAmp)*2.*Pi*PlungeOmega*cos(2*Pi*PlungeOmega*time)

 XRdot(1,1)=-(XfiAmp)*4.*(Pi**2)*(PlungeOmega**2)*sin(2*Pi*PlungeOmega*time)
 XRdot(2,1)=-(YfiAmp)*4.*(Pi**2)*(PlungeOmega**2)*sin(2*Pi*PlungeOmega*time)
 XRdot(3,1)=-(ZfiAmp)*4.*(Pi**2)*(PlungeOmega**2)*sin(2*Pi*PlungeOmega*time)

 elseif(PlugType.eq.2)then
 Xfi = (XfiAmp)*cos(2*Pi*PlungeOmega*time)
 Yfi = (YfiAmp)*cos(2*Pi*PlungeOmega*time)
 Zfi = (ZfiAmp)*cos(2*Pi*PlungeOmega*time)

 MemOri(1,1) = Xfi
 MemOri(2,1) = Yfi
 MemOri(3,1) = Zfi

 Xfidot = -(XfiAmp)*2*Pi*PlungeOmega*sin(2*Pi*PlungeOmega*time)
 Yfidot = -(YfiAmp)*2*Pi*PlungeOmega*sin(2*Pi*PlungeOmega*time)
 Zfidot = -(ZfiAmp)*2*Pi*PlungeOmega*sin(2*Pi*PlungeOmega*time)

 XRdot(1,1)= -(XfiAmp)*4*(Pi**2)*(PlungeOmega**2)*cos(2*Pi*PlungeOmega*time)
 XRdot(2,1)= -(YfiAmp)*4*(Pi**2)*(PlungeOmega**2)*cos(2*Pi*PlungeOmega*time)
 XRdot(3,1)= -(ZfiAmp)*4*(Pi**2)*(PlungeOmega**2)*cos(2*Pi*PlungeOmega*time)
 else
    print *, "Error in reading PlugType"
    stop
 endif
   
 endif

   if(LociSTREAMRank==0.and.LociSTREAMitfsi==1) write(340,'(1x,9(e14.8,1x))') RotVec(1:3),RotVecdot(1:3),RotVecddot(1:3)
   if(LociSTREAMRank==0.and.LociSTREAMitfsi==1) write(345,'(1x,9(e14.8,1x))') MemOri(1:3,1),Xfidot,Yfidot,Zfidot,XRdot(1:3,1)

 end subroutine nlams_undeformed_step
end module nlams_undeformed
