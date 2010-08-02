program main
    implicit none
    integer:: i,j,k
    integer:: id
    integer::nx,ny
    
    real*8    
    real(8)::dx,dy
    real(8)::point(100,2)

    dx=wingspan/(nx-1)
    dy=wingchord/(ny-1)

    
    do j=1,ny
     do i=1,nx
      id = (j-1)*nx+i
      point(id,1) = (i-1)*dx
      point(id,2) = (j-1)*dy
     enddo
    enddo

! connectivity
    id = 0
    do j=1,ny-1
     do i=1,nx-1
        id = (j-1)*ny + i
        elem = 2*id
        if(mod(id,2).ne.0)then
         connect(elem,1) = i
         connect(elem,2) = i+1
         connect(elem,3) = (j-1)*nx + i
        else
         connect(elem,1) = i+1
         connect(elem,2) = (j-1)*nx + i + 1
         connect(elem,3) = (j-1)*nx + i
        endif
      enddo
     enddo
   

!   open(10,file='testrimesh.dat')
!   rewind(10)
!   write(10,*) 'VARIABLES = "X", "Y"'
!   write(10,*) 'ZONE F=FEPOINT, ET=TRIANGLE, N=',NumNodes,',E=',NumElems  
    
