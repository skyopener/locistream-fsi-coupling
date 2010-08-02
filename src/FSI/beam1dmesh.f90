module mesh
    IMPLICIT NONE

		integer :: rank
		integer :: nNodes
		integer :: nElements
		integer :: dofFree
		integer :: nBC
		integer :: dofNode
		integer :: nodePerElement
		integer :: dofElement
    real, dimension(:), allocatable :: dx ! length of each element
    real, dimension(:), allocatable :: nodes 
    real, dimension(:), allocatable :: EI
    real, dimension(:), allocatable :: rhoA
    integer, dimension(:,:), allocatable :: nodeToEquation ! reference from dof nodes to equation nr
    integer, dimension(:,:), allocatable :: elements

    contains
        subroutine assemble_mesh(EIdefault, rhoAdefault)
    
            ! parameters

            ! passed variables
          !  integer :: nNodes, nElements, nodePerElement, dofNode
            real :: EIdefault, rhoAdefault
      !      integer, dimension(:,:) :: bcIndex
            
            ! local variables
            integer :: n, e, i, counter
      !      integer :: nBC

            !nBC = 2 !size(bcIndex, 2)

            ! allocate
            !allocate(nodes(nNodes)) 
            allocate(elements(nodePerElement, nElements), EI(nElements), rhoA(nElements), nodeToEquation(dofNode, nNodes), dx(nElements))
            elements=0; nodeToEquation=0; EI=0.0d0; rhoA=0.0d0; dx=0.0d0
						if (rank==0) print*, '[I] EulerBeam1d: mesh -> mesh succesfully allocated, nElements ', nElements
						!if (rank==0) print*, '[I] nodes', nodes
						
            ! element stiffness matrix
            do e=1,nElements
            	dx(e) = nodes(e+1)-nodes(e)
            end do
            if (rank==0) print*, '[I] EulerBeam1d: mesh -> dx succesfully assembled'

            ! generate mesh
!             nodes = xStart
!             do n=2,nNodes
!                 nodes(n) = nodes(n) + (n-1) * dx(n-1)
!             end do

            ! assemble element information
            do e=1,nElements
                do i=1,nodePerElement
                    elements(i,e) = e + i - 1 ! for a element with two nodes
                end do
            end do
            if (rank==0) print*, '[I] EulerBeam1d: mesh -> elements succesfully assembled'

            ! material properties
            EI = EIdefault
            rhoA = rhoAdefault

            ! equation numbering
            nodeToEquation(1,1) = -1 ! displacement given at x=xStart
            nodeToEquation(2,1) = -1 ! slope given at x=xStart
!             do i=1,nBC
!                 nodeToEquation(bcIndex(1,i), bcIndex(2,i)) = -1
!             end do
						if (rank==0) print*, '[I] EulerBeam1d: mesh -> nodeToEquation1 succesfully assembled'
						
            ! assemble node to equation numbers
    !        if (rank==0) print*, '[I] EulerBeam1d: mesh -> nodeToEquation=', nodeToEquation
            counter = 1
            do n=1,nNodes
                do i=1,dofNode
                    if (nodeToEquation(i,n) /= -1) then
                        nodeToEquation(i,n) = counter
                        counter = counter + 1
                    end if
                end do
            end do
						if (rank==0) print*, '[I] EulerBeam1d: mesh -> nodeToEquation2 succesfully assembled'
!						if (rank==0) print*, '[I] EulerBeam1d: mesh -> nodeToEquation=', nodeToEquation
						
        end subroutine assemble_mesh
        subroutine deallocate_mesh
            
            deallocate(dx, nodeToEquation, elements, rhoA, EI, nodes)            

        end subroutine deallocate_mesh


end module mesh
