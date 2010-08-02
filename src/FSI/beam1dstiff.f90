module stiff
    IMPLICIT NONE
    contains
        subroutine assemble_matrices(sol, q, Pn, M, K, F)
        use mesh
    
        ! parameters

        ! passed Arrays
        real, dimension(:,:) :: sol
        
        real :: q, dw
        real, dimension(:) :: Pn

        real, dimension(:,:) :: M
        real, dimension(:,:) :: K
        real, dimension(:) :: F

        ! local variables
        real, dimension(:,:), allocatable :: ke, ke_base ! element and global stiffness matrix
        real, dimension(:), allocatable :: qe, qe_base  ! element distributed and point force vector and global force vector
        real, dimension(:,:), allocatable :: Me, Me_base  ! element distributed and point force vector and global force vector
        !integer :: dofNode, nElements, nodePerElement, dofElement, nNode
        integer :: e, ni, ii, nj, jj, gi, gj, li, lj

        !nNode = size(Pn, 1)
        !nElements = size(elements, 2)
        !nodePerElement = size(elements, 1)
        !dofNode = size(nodeToEquation, 1)
        !dofElement = nodePerElement * nNode

        allocate(me(dofElement, dofElement), me_base(dofElement, dofElement), ke(dofElement, dofElement), ke_base(dofElement, dofElement), qe(dofElement), qe_base(dofElement))
        ke = 0.; qe=0.; me=0.

        M=0.; K=0.; F=0.

        ke_base = reshape( (/ 12.0d0, -6.0d0, -12.0d0, -6.0d0, -6.0d0, 4.0d0, 6.0d0, 2.0d0, -12.0d0, 6.0d0, 12.0d0, 6.0d0, -6.0d0, 2.0d0, 6.0d0, 4.0d0 /), (/4,4/))
        me_base = reshape( (/ 156.0d0, -22.0d0, 54.0d0, 13.0d0, -22.0d0, 4.0d0, -13.0d0, -3.0d0, 54.0d0, -13.0d0, 156.0d0, 22.0d0, 13.0d0, -3.0d0, 22.0d0, 4.0d0 /), (/4,4/))
        qe_base = (/ 0.5d0, -1.0d0/12.0d0 , 0.50d0, 1.0d0/12.0d0 /)

        ! loop through each element and assemble
        do e=1,nElements
            ! for each element
            me = me_base
            me(2,1) = me(2,1) * dx(e) ;  me(4,1) = me(4,1) * dx(e) ;  me(1,2) = me(1,2) * dx(e) ; me(2,2) = me(2,2) * dx(e)**2 ;  me(3,2) = me(3,2) * dx(e)
            me(4,2) = me(4,2) * dx(e)**2 ; me(2,3) = me(2,3) * dx(e) ;  me(4,3) = me(4,3) * dx(e) ;  
            me(1,4) = me(1,4) * dx(e) ;  me(2,4) = me(2,4) * dx(e)**2 ; me(3,4) = me(3,4) * dx(e) ;  me(4,4) = me(4,4) * dx(e)**2 ;  me = me * rhoA(e) * dx(e) / 420.0d0
            ke = ke_base
            ke(2,1) = ke(2,1) * dx(e) ;  ke(4,1) = ke(4,1) * dx(e) ;  ke(1,2) = ke(1,2) * dx(e) ; ke(2,2) = ke(2,2) * dx(e)**2 ;  ke(3,2) = ke(3,2) * dx(e)
            ke(4,2) = ke(4,2) * dx(e)**2 ; ke(2,3) = ke(2,3) * dx(e) ;  ke(4,3) = ke(4,3) * dx(e) ;  
            ke(1,4) = ke(1,4) * dx(e) ;  ke(2,4) = ke(2,4) * dx(e)**2 ; ke(3,4) = ke(3,4) * dx(e) ;  ke(4,4) = ke(4,4) * dx(e)**2 ;  ke = ke * EI(e) / dx(e)**3
            qe = qe_base
	    qe = qe * dx(e)
            qe(2) = qe(2)*dx(e) ;  qe(4) = qe(4)*dx(e)
            !qe = q * dx(e) * qe
	    qe = 0.5d0 * ( Pn(elements(1,e)) + Pn(elements(2,e) ) ) * qe ! Pn / dx * dx * qe
	    dw = sol(1, elements(2,e)) - sol(1, elements(1,e)) ! difference in w at node 2 - node 1
	    qe = qe * ( dx(e) / sqrt( dx(e)**2 + dw**2 ) ) ! taking the y component
            do ni=1,nodePerElement
            do ii=1,dofNode
                if ( (nodeToEquation(ii,elements(ni,e)) /= -1) ) then
                    !F( nodeToEquation(ii,elements(ni,e))) = F( nodeToEquation(ii,elements(ni,e)) ) + qe( dofNode*(ni-1)+ii ) + Pn( elements(ni,e) ) * (2.0d0 - real(ii)) ! Pn only if ii=1
		    F( nodeToEquation(ii,elements(ni,e))) = F( nodeToEquation(ii,elements(ni,e)) ) + qe( dofNode*(ni-1)+ii ) ! Pn represented as distributed loading q
                end if
                do nj=1,nodePerElement
                do jj=1,dofNode
                    if ( (nodeToEquation(ii,elements(ni, e)) /= -1) .AND. (nodeToEquation(jj,elements(nj,e)) /= -1) ) then
                        gi = nodeToEquation(ii,elements(ni,e)); gj = nodeToEquation(jj,elements(nj,e))
                        li = dofNode*(ni-1)+ii ; lj = dofNode*(nj-1)+jj
        !                print*, gi, gj, e, li, lj
                        K( gi, gj ) = K( gi, gj ) + ke(li, lj)
                        M( gi, gj ) = M( gi, gj ) + me(li, lj)
                    end if
                    if ( (nodeToEquation(ii,elements(ni,e)) /= -1) .AND. (nodeToEquation(jj,elements(nj,e)) == -1) ) then
                        gi = nodeToEquation(ii,elements(ni,e)); gj = nodeToEquation(jj,elements(nj,e))
                        li = dofNode*(ni-1)+ii ; lj = dofNode*(nj-1)+jj
                        F(gi) = F(gi) - ke(li, lj) * sol(jj, nj)
!                        F(gi) = F(gi) - ke(li, lj) * sol(jj, nj)
                    end if
                end do
                end do
            end do
            end do
        end do

        deallocate( ke, ke_base, qe_base, qe, me, me_base)

        end subroutine assemble_matrices
end module stiff
