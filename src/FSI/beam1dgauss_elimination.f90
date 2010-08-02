module gauss
    IMPLICIT NONE
    contains
        subroutine lin_solve(A, b)
        use mesh
    
        ! parameters
        real, parameter :: ONE = 1.0d0
        real, parameter :: SMALL = 1.0d-10

        ! passed Arrays
        real, dimension(:,:) :: A
        real, dimension(:) :: b

        ! local variables
        integer :: N
        integer :: i, k
        real :: m

        ! initialization
        N = dofFree !size(b)

        !------------------------------------------------------------------------------
        ! inverse start
        ! Reduce augmented matrix to upper trianguar matrix
        ! Remember for this project, diagonal always 1.0
        do k=1, N-1 ! sweep through rows
            do i=k+1, N
            		!if (rank.eq.0) print*, 'N = ', N, 'k = ', k, 'i = ', i, 'A(k,i) = ', A(k,i)
                if (abs(A(k,i)) > SMALL ) then
                		!if (rank.eq.0) print*, 'N = ', N, 'k = ', k, 'A(k,k) = ', A(k,k), 'A(k,i) = ', A(k,i)
                    m = A(k,i) / (A(k,k) )
                    !if (rank.eq.0) print*, 'After m = ', m
                    A(:,i) = A(:,i) - m * A(:, k)
                    !if (rank.eq.0) print*, 'After A(:,i) = '
                    b(i) = b(i) - m * b(k)
                    !if (rank.eq.0) print*, 'After b(i) = '
                end if
            end do
        end do

				!if (rank.eq.0) print*, 'Before diagonal'
        ! Make diagonals one
        do i=1,N
            m = A(i,i) 
            A(:,i) = A(:,i) / m
            b(i) = b(i) / m
        end do
        !if (rank.eq.0) print*, 'after diagonal'

        ! Backward sweep
        do k=N, 2, -1 ! sweep through rows
            do i=k-1, 1, -1
                if (abs(A(k,i)) > SMALL ) then
!                    m = A(i,k) / A(k,k)
!                    A(i,:) = A(i,:) - m * A(k, :)
                    b(i) = b(i) - A(k,i) * b(k)
                end if
            end do
        end do
				!if (rank.eq.0) print*, 'after backsub'
				
        end subroutine lin_solve
end module gauss
