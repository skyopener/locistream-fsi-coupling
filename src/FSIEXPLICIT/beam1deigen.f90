module eigen
    IMPLICIT NONE
    contains
        subroutine compute_eigen_lapack(a, d, v)
        ! passed variables
        real, dimension(:), intent(inout) :: d
        real, dimension(:,:), intent(inout) :: a
        real, dimension(:,:), intent(inout) :: v
        
        ! local
        integer :: n, lda, ldvr, lwmax, lwork, info
        real, dimension(:,:), allocatable :: dummy
        real, dimension(:), allocatable :: work, dummyv
        n = size(a,1)
        lda = n
        ldvr = size(v,1)
        lwmax = 100000

        allocate(work(lwmax), dummy(1,1), dummyv(n))

        ! query the optimal workspace
        lwork = -1
!        call dgeev('n', 'v', n, a, lda, d, dummyv, dummy, 1, v, ldvr, work, lwork, info)
!        lwork = min( lwmax, int(work(1)))
!        print*, 'work 1 = ', work(1)
!        print*, 'lwork = ', lwork
        lwork = lwmax

!        print*, 'a = ', a

        !solve
        call dgeev('n', 'v', n, a, lda, d, dummyv, dummy, 1, v, ldvr, work, lwork, info)

!        print*, 'a = ', a
!        print*, 'd = ', d
!        print*, 'di = ', dummyv
!        print*, 'v = ', v

        ! check
        if (info>0) then
            print*, 'The algorithm failed to computed eigenvalues.'
        end if

        deallocate(work, dummy, dummyv)

        end subroutine compute_eigen_lapack



        subroutine compute_eigen(a, d, v, nrot)
        ! Computes all eigenvalues and eigenvectors of a real symmetric NxN
        ! matrix a. On output, elements of a above the diagonal are destroyed.
        ! d : is a vector of length N that returns the eigenvalues of a
        ! v : an N x N matrix whose columns contain, on output, the normalized eigen vectors of a
        ! nrot : returns the number of Jacobi rotations that were required
    
        ! passed variables
        integer, intent(out) :: nrot
        real, dimension(:), intent(out) :: d
        real, dimension(:,:), intent(inout) :: a
        real, dimension(:,:), intent(out) :: v

        ! local variables
        integer :: i,j, ip, iq, n
        real :: c, g, h, s, sm, t, tau, theta,tresh
        real, dimension(size(d)) :: b, z

        ! get n, if not all the same report error
!        n = size(a,1)
        if ((size(a,1) == size(a,2)) .and.  (size(v,1) == size(v,2)) .and. (size(a,1) == size(v,1)) .and. (size(a,1) == size(d) ) ) then
            n = size(a,1)
        else
            print*, 'size error'
        end if
        print*, 'n = ', n

        ! initialize v to the identity matrix
        v = 0.0
        do i=1,n
            v(i,i) = 1.d0
        end do
        ! initialize b and d to be the diagonal of a
        do i=1,n
            b(i) = a(i,i)
        end do
        d = b
        z = 0.0
        nrot = 0
        do i=1,1000
                print*, 'nrot = ', nrot
            ! sum of off-diagonal elements
            sm = 0.0
            do j=1,n-1
                sm = sm + sum(abs(a(j,j+1:n))) 
            end do
            if (sm == 0.0) then
                print*, 'nrot = ', nrot
                return ! the normal return, which relies on quadratic convergence to machine underflow
            end if
            tresh = merge(0.2 * sm / n**2, 0.0, i<4) ! on the first three sweeps, we will rotate only if tresh exceeded
            do ip=1, n-1
                do iq=ip+1,n
                    g = 100.0 * abs(a(ip,iq)) 
                    ! after 4 sweeps, skip the rotation if the off-diagonal element is small
                    if ((i>4) .and. (abs(d(ip))+g == abs(d(ip))) .and. (abs(d(iq))+g == abs(d(ip)))) then
                        a(ip,iq) = 0.0
                    else if (abs(a(ip,iq)) > tresh ) then
                        h = d(iq) - d(ip)
                        if (abs(h)+g == abs(h)) then
                            t = a(ip,iq) / h
                        else
                            theta = 0.5 * h / a(ip,iq)
                            t = 1.0 / (abs(theta) + sqrt(1.0+theta**2))
                            if (theta < 0.0) t = -t
                        end if

                        c = 1.0 /sqrt(1.+t**2)
                        s = t*c
                        tau=s/(1.0 + c)
                        h = t*a(ip,iq)
                        z(ip)=z(ip)-h
                        z(iq)=z(iq)+h
                        d(ip)=d(ip)-h
                        d(iq)=d(iq)+h
                        a(ip,iq)=0.0
                        call jrotate(a(1:ip-1,ip),a(1:ip-1,iq)) ! case of rotations 1<= j < p
                        call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq)) ! case of rotations p<= j < q
                        call jrotate(a(ip,iq+1:n),a(iq,iq+1:n)) ! case of rotations q<j <= n
                        call jrotate(v(:,ip),v(:,iq)) 
                        nrot=nrot+1
                    end if
                end do
            end do
            b=b+z
            d=b
            z=0.0
         end do
         print*, 'iterations over 100 needed!!'

         contains

         subroutine jrotate(a1,a2)
            real, dimension(:), intent(inout) :: a1, a2
            real, dimension(size(a1)) :: wk1

            wk1 = a1
            a1 = a1 - s*(a2 + a1 * tau)
            a2 = a2 + s*(wk1 - a2 * tau)
        end subroutine jrotate
    end subroutine compute_eigen
    subroutine sort_eigen(d, v)

        real, dimension(:), intent(inout) :: d
        real, dimension(:,:), intent(inout) :: v

        integer :: i, j, n
        real :: dummy
        real, dimension(size(d)) :: dummy_v

        n = size(d)
        do i=1,n-1
            j = imaxloc(d(i:n)) + i - 1
            if (j /= i) then
                !call swap(d(i), d(j))
                dummy = d(i)
                d(i) = d(j)
                d(j) = dummy
                !call swap(v(:,i), v(:,j))
                dummy_v = v(:,i)
                v(:,i) = v(:,j)
                v(:,j) = dummy_v
            end if
        end do

        contains
        function imaxloc(arr)
            real, dimension(:), intent(in) :: arr
            integer :: imaxloc
            integer, dimension(1) :: imax
            imax = maxloc(arr(:))
            imaxloc = imax(1)
        end function imaxloc
    end subroutine sort_eigen
    
end module eigen
