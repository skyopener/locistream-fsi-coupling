module MATRIX_INVERTER

    IMPLICIT NONE
    contains

    function inverse(a)

        ! passed variables 
        real, dimension(:,:), intent(inout) :: a

        ! output
        real, dimension(:,:), allocatable :: inverse


        !	Local Arrays
        real ::   m, SMALL
        INTEGER ::  n, i, j, k


        !	Initialize parameters
        n = size(a,1)
        SMALL = 1.D-12

        !       Allocate matrices
        ALLOCATE( inverse(n, n) )
        inverse=0.

        ! Fill in the augmented part with identity matrix
        DO i=1, n
                inverse(i, i) = 1.d0
        ENDDO

        !------------------------------------------------------------------------------
        ! inverse start
        ! Reduce augmented matrix to upper trianguar matrix
        ! Remember for this project, diagonal always 1.0
        !IF (1 == 0) THEN
        DO k=1, n-1 ! sweep through rows
                DO i=k+1,n
                        IF (abs(A(i,k)) > SMALL) THEN
                                m = A(i,k) / A(k,k)
                                A(i,:) = A(i,:) - m * A(k,:)
                                inverse(i,:) = inverse(i,:) - m * inverse(k,:)
                        END IF
                ENDDO
        END DO

        ! Make diagonals one
        DO i=1,n
                m = A(i,i)
                A(i,:) = A(i,:) / m
                inverse(i,:) = inverse(i,:) / m
        END DO

        ! Backward sweep
        DO k=n, 2, -1 ! sweep through rows
                DO i=k-1,1,-1
                        IF (abs(A(i,k)) > SMALL) THEN
                                m = A(i,k) / A(k,k) ! PIVOTROW(k) should be one!
                                A(i,:) = A(i,:) - m * A(k,:)
                                inverse(i,:) = inverse(i,:) - m * inverse(k,:)
                        END IF
                ENDDO
        END DO

    END function

end module MATRIX_INVERTER


