!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Module.- NLAMS_FORM_SHELL
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!     
!  
!-> Reference:
!     Paper: Corotational Nonlinear Analysis of Thin Plates and Shells Using a New Shell Element 
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 69, pp. 859-885, 2007
!     The above reference refers to another journal paper whose details are below:
!     Paper: A Study of Three-Node Triangular Plate Bending Elements
!     Journal: International Journal of Numerical Methods in Engineering, Vol. 15, pp. 1771-1812, 1980.
!
!->Subroutines:
!
!     |-nlams_shell
!
!
!->Remarks:
!    Variables are precisely as presented in the reference paper.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nlams_form_shell

contains
subroutine nlams_shell(Bm,Bb,K_elem_opt,K_elem_dkt,K_elem_shell)
!subroutine nlams_shell(Bb,K_elem_opt,K_elem_dkt,K_elem_shell)
use nlams_comp_abbd
implicit none
!intent(in) variables
real(8),dimension(3,9),intent(in)::Bb,Bm
real(8),dimension(9,9),intent(in)::K_elem_opt,K_elem_dkt

!intent(out) variables
real(8),dimension(18,18),intent(out)::K_elem_shell

!Local variables
integer :: i,j,p,q
real(8),dimension(9,9)::K_elem_shell_12,K_elem_shell_21

K_elem_shell_12 = 0.0
K_elem_shell_21 = 0.0

! Eq. (33) in referece paper 
! K_elem_shell_12: correspoding to Kmb in Eq. (33)
! K_elem_shell_21: correspoding to kbm in Eq. (33)

K_elem_shell_12 =  K_elem_shell_12 + matmul(matmul(transpose(Bm),Be_ij),Bb)
K_elem_shell_21 =  K_elem_shell_21 + matmul(matmul(transpose(Bb),Be_ij),Bm)

! if Be_ij is zero, K_elem_shell_12, K_elem_shell_21 = 0.0
!write(*,*) 'kmb'
!do i =1, 9
!   write(*,'(9e14.6)') (K_elem_shell_12(i,j),j=1,9)
!enddo
! write(*,*) 'kbm'
!do i =1, 9
!   write(*,'(9e14.6)') (K_elem_shell_21(i,j),j=1,9)
!enddo
!read(*,*)

! K_elem_shell(1:9,1:9) = Km
do i=1,9
 do j=1,9
      K_elem_shell(i,j) = K_elem_opt(i,j)
 enddo
enddo

! i=1,9, j=9,18: kmb
p = 0; q = 0
do i=1,9
   p = p + 1
   do j=10,18
      q = q + 1
      K_elem_shell(i,j) = K_elem_shell_12(p,q)
   enddo
   q = 0
enddo

! i=10,18, j=1,9: kbm
p = 0;q = 0
do i=10,18
   p = p + 1
   do j=1,9
      q = q + 1
      K_elem_shell(i,j) = K_elem_shell_21(p,q)
   enddo
   q = 0
enddo

! i=10,18, j=10,18: Kb
p = 0;q = 0
do i=10,18
   p = p + 1
   do j=10,18
      q = q + 1
      K_elem_shell(i,j) = K_elem_dkt(p,q)
   enddo
   q = 0
enddo

!Check stiffness matrix of the shell
!write(*,*) 'Stiffness matrix of shell[18*18] '
!do i =1, 18
!   write(*,'(18e12.4)') (K_elem_shell(i,j),j=1,18)
!enddo
!read(*,*)

 end subroutine nlams_shell
end module nlams_form_shell


