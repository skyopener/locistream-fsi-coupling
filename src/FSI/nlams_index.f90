!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-> Copyright by The University of Michigan, Aerospace Department. 2008
!
!-> Subroutine.- NLAMS_index. Satish Chimakurthi.
!
!-> Language: FORTRAN90, Free Format.
!
!-> Description:
!
!   Generates the DOF iConnectices for all the nodes of each element
!
!-> Reference:
!
!
!->Subroutines:
!
!     |- nlams_index
!
!->Remarks:
!  1) NumNodesPerEl is number of nodes per element
!  2) el is element number whose system dofs are to be determined
!  3) NumDofPerNode is number of dofs per node^M
!  4) Connect is the connectivity array 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine nlams_index(iElem)
 use nlams_init, only: ElDofIndex
 use mesh_connect_vars,only:NumNodesPerEl,NumDofPerNode,Connect

 implicit none

 integer::edof,i,j,start,tt
 integer,dimension(3)::nd
 integer,intent(in)::iElem

 edof = NumNodesPerEl*NumDofPerNode
 
 nd = Connect(iElem,1:3)

 tt=0
 do i=1,NumNodesPerEl
    start=(nd(i)-1)*NumDofPerNode
    do j=1,NumDofPerNode
        tt= tt + 1
        ElDofIndex(tt)= start + j
    enddo
 enddo

end subroutine nlams_index

