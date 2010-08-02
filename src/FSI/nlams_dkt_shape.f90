!Shape functions and their derivatives for the DKT

 subroutine dkt_shape (Gauss)
   use ElemLocal
   use nlams_init, only:N_dkt,dN_dkt_g2,dN_dkt_g3

   implicit none

   real(8),dimension(3),intent(in)::Gauss

!The shape functions below are exactly as in Eq. 39 of the reference paper of
!Khosravi et al. mentioned above
    N_dkt(1,1) = ((Gauss(1)**2)*(3.d0 - 2.d0*Gauss(1))) + (2.d0*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,2) = (-Gauss(1)**2)*(y12*Gauss(2) + y13*Gauss(3)) - ((0.5d0)*(y12 + y13)*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,3) = ( Gauss(1)**2)*(x12*Gauss(2) + x13*Gauss(3)) + ((0.5d0)*(x12 + x13)*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,4) = ((Gauss(2)**2)*(3.d0 - 2.d0*Gauss(2))) + (2.d0*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,5) = (-Gauss(2)**2)*(y23*Gauss(3) + y21*Gauss(1)) - ((0.5d0)*(y23 + y21)*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,6) = ( Gauss(2)**2)*(x23*Gauss(3) + x21*Gauss(1)) + ((0.5d0)*(x23 + x21)*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,7) = ((Gauss(3)**2)*(3.d0 - 2.d0*Gauss(3))) + (2.d0*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,8) = (-Gauss(3)**2)*(y31*Gauss(1) + y32*Gauss(2)) - ((0.5d0)*(y31 + y32)*Gauss(1)*Gauss(2)*Gauss(3))
    N_dkt(1,9) = ( Gauss(3)**2)*(x31*Gauss(1) + x32*Gauss(2)) + ((0.5d0)*(x31 + x32)*Gauss(1)*Gauss(2)*Gauss(3))
    
    !Derivatives of shape functions w.r.t. Gauss(2) and Gauss(3)
    
    !w.r.t. Gauss(2)
    dN_dkt_g2(1,1) =  (2.d0*(1.d0-Gauss(2)-Gauss(3))**2)  - (2.d0*(1.d0+2.d0*Gauss(2)+2.d0*Gauss(3)))*(1.d0-Gauss(2)-Gauss(3)) - (2.d0*Gauss(2)*Gauss(3)) + (2.d0*Gauss(3)*(1.d0 - Gauss(2) - Gauss(3)))
    dN_dkt_g2(1,2) = -(2.d0*Gauss(2) - 2.d0 + 2.d0*Gauss(3))*(y12*Gauss(2) + y13*Gauss(3)) - y12*(1.d0 - Gauss(2)-Gauss(3))**2 - 0.5d0*(y12+y13)*(1.d0-Gauss(2)-Gauss(3))*Gauss(3) + 0.5d0*Gauss(2)*Gauss(3)*(y12+y13)
    dN_dkt_g2(1,3) =  (2.d0*Gauss(2) - 2.d0 + 2.d0*Gauss(3))*(x12*Gauss(2) + x13*Gauss(3)) + x12*(1.d0 + Gauss(2)**2 - 2.d0*Gauss(2) + Gauss(3)**2 - 2.d0*Gauss(3) + 2.d0*Gauss(2)*Gauss(3)) + 0.5d0*(x12+x13)*Gauss(3)*(1.d0-Gauss(2)-Gauss(3)) - 0.5*(x12+x13)*Gauss(2)*Gauss(3)
    dN_dkt_g2(1,4) =  2.d0*Gauss(2)*(3.d0-2.d0*Gauss(2)) - 2.d0*Gauss(2)**2 + 2.d0*Gauss(3)*(1.d0-Gauss(2)-Gauss(3)) - 2.d0*Gauss(2)*Gauss(3)
    dN_dkt_g2(1,5) = -2.d0*Gauss(2)*(y23*Gauss(3) + y21 - y21*Gauss(2) - y21*Gauss(3)) + y21*Gauss(2)**2 - 0.5d0*(y23+y21)*(1.d0-Gauss(2)-Gauss(3))*Gauss(3) + 0.5d0*Gauss(2)*Gauss(3)*(y23+y21) 
    dN_dkt_g2(1,6) =  2.d0*Gauss(2)*(x23*Gauss(3) + x21 - x21*Gauss(2) - x21*Gauss(3)) - (x21*Gauss(2)**2) + 0.5d0*(x23+x21)*Gauss(3)*(1.d0-Gauss(2)-Gauss(3)) - 0.5d0*Gauss(2)*Gauss(3)*(x23+x21)
    dN_dkt_g2(1,7) =  2.d0*Gauss(3)*(1.d0-Gauss(2)-Gauss(3))-(2.d0*Gauss(2)*Gauss(3))
    dN_dkt_g2(1,8) =  (-Gauss(3)**2)*(-y31+y32) - 0.5d0*(y31+y32)*(1.d0-Gauss(2)-Gauss(3))*Gauss(3) + 0.5d0*(y31+y32)*Gauss(2)*Gauss(3)
    dN_dkt_g2(1,9) =  ( Gauss(3)**2)*(-x31+x32) + 0.5d0*(x31+x32)*(1.d0-Gauss(2)-Gauss(3))*Gauss(3) - 0.5d0*(x31+x32)*Gauss(2)*Gauss(3)
    
    !w.r.t. Gauss(3)
    dN_dkt_g3(1,1) =  2.d0*(1.d0+Gauss(2)**2 - 2.d0*Gauss(2) + Gauss(3)**2 - 2.d0*Gauss(3)*(1.d0-Gauss(2))) + (1.d0+2.d0*Gauss(2)+2.d0*Gauss(3))*(2.d0*Gauss(3) - 2.d0 + 2.d0*Gauss(2)) + (2.d0*Gauss(2)*(1.d0-Gauss(2)-Gauss(3))) - 2.d0*Gauss(2)*Gauss(3)
    dN_dkt_g3(1,2) = -(2.d0*Gauss(3) - 2.d0*(1-Gauss(2)))*(y12*Gauss(2) + y13*Gauss(3)) - (1.d0+Gauss(2)**2 - 2.d0*Gauss(2) + Gauss(3)**2 - 2.d0*(1.d0-Gauss(2))*Gauss(3))*y13 - 0.5d0*(y12+y13)*Gauss(2)*(1.d0-Gauss(2)-Gauss(3)) + 0.5d0*(y12+y13)*Gauss(2)*Gauss(3)
    dN_dkt_g3(1,3) =  x13*(1.d0+Gauss(2)**2 - 2.d0*Gauss(2) + Gauss(3)**2 - 2.d0*Gauss(3) + 2.d0*Gauss(2)*Gauss(3)) + (x12*Gauss(2)+x13*Gauss(3))*(2.d0*Gauss(3)-2.d0+2.d0*Gauss(2)) + 0.5d0*(x12+x13)*Gauss(2)*(1.d0-Gauss(2)-Gauss(3)) - 0.5d0*(x12+x13)*Gauss(2)*Gauss(3)
    dN_dkt_g3(1,4) =  2.d0*(1.d0-Gauss(2)-Gauss(3))*Gauss(2) - 2*Gauss(2)*Gauss(3)
    dN_dkt_g3(1,5) =  (-(y23-y21)*Gauss(2)**2)  - 0.5d0*(y23+y21)*(1.d0-Gauss(2)-Gauss(3))*Gauss(2) + (0.5d0*(y23+y21)*Gauss(2)*Gauss(3))
    dN_dkt_g3(1,6) =  (Gauss(2)**2)*(x23-x21)   + 0.5d0*(x23+x21)*Gauss(2)*(1.d0-Gauss(2)-Gauss(3)) - (0.5d0*Gauss(2)*Gauss(3)*(x23+x21))
    dN_dkt_g3(1,7) =  2.d0*Gauss(3)*(3.d0-2.d0*Gauss(3)) - 2.d0*Gauss(3)**2 + 2.d0*Gauss(2)*(1.d0-Gauss(2)-Gauss(3)) - 2.d0*Gauss(2)*Gauss(3)
    dN_dkt_g3(1,8) = -2.d0*Gauss(3)*(y31-y31*Gauss(2)-y31*Gauss(3) + y32*Gauss(2)) + (y31*Gauss(3)**2) - 0.5d0*(y31+y32)*Gauss(2)*(1.d0-Gauss(2)-Gauss(3)) + (0.5d0*(y31+y32)*Gauss(2)*Gauss(3))
    dN_dkt_g3(1,9) =  2.d0*Gauss(3)*(x31-x31*Gauss(2)-x31*Gauss(3) + x32*Gauss(2)) - (x31*Gauss(3)**2) + 0.5d0*(x31+x32)*Gauss(2)*(1.d0-Gauss(2)-Gauss(3)) - (0.5d0*(x31+x32)*Gauss(2)*Gauss(3))
  
  end subroutine dkt_shape

