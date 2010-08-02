!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Subroutine.- opt_shape. Satish Chimakurthi.
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   Compute shape functions and their derivative of OPT element
!
!-> Reference:
!  Corotational Nonlinear Analysis of Thin Plates and Shells Using a
!  New Shell Element, Journal: International Journal of Numerical Methods in
!  Engineering, Vol.69, pp. 859-885, 2007
!
!->Subroutines:
!
!     |- opt_shape
!
!->Remarks:
!  1) NumNodesPerEl is number of nodes per element
!  2) el is element number whose system dofs are to be determined
!  3) NumDofPerNode is number of dofs per node^M
!  4) Connect is the connectivity array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine opt_shape (Gauss)
 use ElemLocal
 use nlams_input, only: alfab
 use nlams_init,  only: N_opt, dN_opt_g2, dN_opt_g3

 implicit none
 real(8),dimension(3),intent(in)::Gauss

!The shape functions below are exactly as in Eq. 40 in the paper of Khosravi et
!al. as indicated above.
 N_opt(1,1) = Gauss(1)
 N_opt(1,2) = 0.0d0
 N_opt(1,3) = 0.5d0*(alfab*Gauss(1))*(y12*Gauss(2) - y31*Gauss(3))  
 N_opt(1,4) = 1.0d0
 N_opt(1,5) = 0.0d0
 N_opt(1,6) = 0.5d0*(alfab*Gauss(2))*(y23*Gauss(3) - y12*Gauss(1))
 N_opt(1,7) = Gauss(3)
 N_opt(1,8) = 0.0d0
 N_opt(1,9) = 0.5d0*(alfab*Gauss(3))*(y31*Gauss(1) - y23*Gauss(2))

 N_opt(2,1) = 0.0d0
 N_opt(2,2) = Gauss(1)
 N_opt(2,3) = 0.5d0*(alfab*Gauss(1))*(x21*Gauss(2) - x13*Gauss(3))
 N_opt(2,4) = 0.0d0
 N_opt(2,5) = Gauss(2)
 N_opt(2,6) = 0.5d0*(alfab*Gauss(2))*(x32*Gauss(3) - x21*Gauss(1))
 N_opt(2,7) = 0.0d0
 N_opt(2,8) = Gauss(3)
 N_opt(2,9) = 0.5d0*(alfab*Gauss(3))*(x13*Gauss(1) - x32*Gauss(2))

!Derivatives of shape functions w.r.t. Gauss(2) and Gauss(3)
!w.r.t Gauss(2)

!u
 dN_opt_g2(1,1) = -1.0d0
 dN_opt_g2(1,2) =  0.0d0
 dN_opt_g2(1,3) = -0.5d0*alfab*(y12*Gauss(2) - y31*Gauss(3)) + 0.5d0*alfab*y12*(1-Gauss(2)-Gauss(3)) 
 dN_opt_g2(1,4) =  1.0d0
 dN_opt_g2(1,5) =  0.0d0
 dN_opt_g2(1,6) =  0.5d0*alfab*Gauss(2)*y12 + 0.5d0*alfab*(y23*Gauss(3) - y12 + Gauss(2)*y12 + Gauss(3)*y12)
 dN_opt_g2(1,7) =  0.0d0
 dN_opt_g2(1,8) =  0.0d0
 dN_opt_g2(1,9) = -0.5d0*alfab*Gauss(3)*(y31 + y23)

!v
 dN_opt_g2(2,1) =  0.0d0
 dN_opt_g2(2,2) = -1.0d0
 dN_opt_g2(2,3) = -0.5d0*alfab*(x21*Gauss(2) - x13*Gauss(3)) + 0.5d0*alfab*x21*(1-Gauss(2)-Gauss(3))
 dN_opt_g2(2,4) =  0.0d0
 dN_opt_g2(2,5) =  1.0d0
 dN_opt_g2(2,6) = (0.5d0*alfab*Gauss(2)*x21) + 0.5d0*alfab*(x32*Gauss(3) - x21 + Gauss(2)*x21 + Gauss(3)*x21)
 dN_opt_g2(2,7) =  0.0d0
 dN_opt_g2(2,8) =  0.0d0
 dN_opt_g2(2,9) = -0.5d0*alfab*Gauss(3)*(x13 + x32)

!w.r.t Gauss(3)
!u
 dN_opt_g3(1,1) = -1.0d0
 dN_opt_g3(1,2) =  0.0d0
 dN_opt_g3(1,3) = -0.5d0*alfab*(y12*Gauss(2) - y31*Gauss(3)) - 0.5d0*alfab*y31*(1-Gauss(2)-Gauss(3))
 dN_opt_g3(1,4) =  0.0d0
 dN_opt_g3(1,5) =  0.0d0
 dN_opt_g3(1,6) =  0.5d0*alfab*Gauss(2)*(y23+y12)
 dN_opt_g3(1,7) =  1.0d0
 dN_opt_g3(1,8) =  0.0d0
 dN_opt_g3(1,9) = (0.5d0*alfab*y31*(1-Gauss(2)-Gauss(3))) - (0.5d0*alfab*y31*Gauss(3)) - 0.5d0*alfab*y23*Gauss(2)

!v
 dN_opt_g3(2,1) =  0.0d0
 dN_opt_g3(2,2) = -1.0d0
 dN_opt_g3(2,3) = -0.50*alfab*(x21*Gauss(2) - x13*Gauss(3)) - 0.5d0*alfab*x13*(1-Gauss(2)-Gauss(3))
 dN_opt_g3(2,4) =  0.0d0
 dN_opt_g3(2,5) =  0.0d0
 dN_opt_g3(2,6) =  0.5d0*(alfab*Gauss(2))*(x32 + x21)
 dN_opt_g3(2,7) =  0.0d0
 dN_opt_g3(2,8) =  1.0d0
 dN_opt_g3(2,9) = -0.5d0*alfab*Gauss(3)*x13 + (0.5d0*alfab*(x13*Gauss(1) - x32*Gauss(2)))

end subroutine opt_shape

