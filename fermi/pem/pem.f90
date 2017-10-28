!Wirte by Liao Yuanda
!####################################################################
module globle_vab
    implicit none
    real*8, parameter :: PI = 3.1415926535898
    real*8 :: a, b, beta, miu, Emax, Emin
end module globle_vab
!********************************************************************
real*8 function logweight_func(x)
    use globle_vab
    implicit none
    real*8 :: x
	!write(*,*) beta, miu
    logweight_func = log( 1.0+exp( -beta*(a*x+b-miu)) )
end function logweight_func
!********************************************************************
real*8 function kernel(cutoff, i)
    use globle_vab
    implicit none
    integer :: cutoff, i

    kernel = 1.0/(cutoff+1.0)*( (cutoff-i+1.0)*cos(PI*i/(cutoff+1.0))+sin(PI*i/(cutoff+1.0))*atan(PI/(cutoff+1.0)) )
end function

!********************************************************************

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!Interface Function
subroutine pem_logweight(cutoff, rank, nonzeros, values, colind, rowptr, logweight)
    use globle_vab
    implicit none
    integer :: cutoff, rank, nonzeros, colind(nonzeros), rowptr(rank+1), i, j
    real*8 :: logweight
    complex*16 :: values(nonzeros)
    real*8, allocatable :: moments(:), coeffs(:)
    real*8, external :: logweight_func

    allocate( moments(cutoff), coeffs(cutoff))
    call pem_calculate_coeffs(cutoff, coeffs, logweight_func)
    !do i = 1, cutoff
	!	write(*,*) coeffs(i)
	!end do
	call pem_calculate_moments(cutoff, rank, nonzeros, values, colind, rowptr, moments)
    call pem_kernel_expansion(cutoff, moments, coeffs, logweight)

    deallocate(moments, coeffs)

end subroutine pem_logweight
!********************************************************************
subroutine pem_calculate_moments(cutoff, rank, nonzeros, values, colind, rowptr, moments)
    implicit none
    integer :: cutoff, rank, nonzeros, colind(nonzeros), rowptr(rank+1), i
    real*8 :: moments(cutoff)
    complex*16 :: values(nonzeros)

    moments = 0.0
    moments(1) = rank
    do i = 1, rank	
        call pem_diagonal_element(cutoff, rank, nonzeros, values, colind, rowptr, moments, i)
    end do

    do i = 3, cutoff-1, 2
        moments(i) = 2.0*moments(i)-moments(1)
    end do
    do i = 4, cutoff-1, 2
        moments(i) = 2.0*moments(i)-moments(2)
    end do
end subroutine pem_calculate_moments
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!********************************************************************
subroutine pem_kernel_expansion(cutoff, moments, coeffs, logweight)
    implicit none
    integer :: i, cutoff
    real*8 :: logweight, moments(cutoff), coeffs(cutoff), kernel

    logweight = 0.0
    do i = 1, cutoff
        logweight = logweight+kernel(cutoff,i)*moments(i)*coeffs(i)
    end do
end subroutine pem_kernel_expansion
!********************************************************************
subroutine pem_calculate_coeffs(cutoff, coeffs, f)
    use globle_vab
    implicit none
    real*8, external :: f
    integer :: cutoff, k, m
    real*8 :: coeffs(cutoff), x, sumover

    do m = 1, cutoff
        sumover = 0.0
        do k = 1, 2*cutoff
            x = cos(PI*(k-0.5)/(2*cutoff))
            sumover = sumover+f(x)*cos(PI*(m-1)*(k-0.5)/(2*cutoff))
        end do
        coeffs(m) = 2.0*sumover/(2*cutoff)
    end do
    coeffs(1) = coeffs(1)/2.0

end subroutine pem_calculate_coeffs
!********************************************************************
subroutine pem_diagonal_element(cutoff, rank, nonzeros, values, colind, rowptr, moments, ket)
    implicit none
    integer :: cutoff, rank, nonzeros, colind(nonzeros), rowptr(rank+1), ket, i, j, m
    real*8 :: moments(cutoff)
    complex*16 :: values(nonzeros), sum1, sum2, keep
    complex*16, allocatable :: tmp(:), jm0(:), jm1(:)

    allocate ( tmp(rank), jm0(rank), jm1(rank) )
    tmp = (0.0,0.0)
    jm0 = (0.0,0.0)
    jm1 = (0.0,0.0)
    jm0(ket) = (1.0,0.0) ! set |j,0>

    ! calculate |j,2> = X|j,1>
    call pem_sparse_product(rank, nonzeros, values, colind, rowptr, jm1, jm0)
    sum1 = conjg(jm0(ket)*jm1(ket))
    moments(2) = moments(2)+real(sum1)

    sum2 = (0.0,0.0)
    do j = 1, rank
        sum2 = sum2+conjg(jm1(j))*jm1(j)
    end do
    moments(3) = moments(3)+real(sum2)

    !calculate |j,m> = 2X|j,m-1> - |j,m-2>
    do m = 2, cutoff/2-1
     ! calculate |tmp> = X|jm1>
        call pem_sparse_product(rank, nonzeros, values, colind, rowptr, tmp, jm1)
        sum1 = (0.0,0.0)
        sum2 = (0.0,0.0)
        do i = 1, rank
            keep = tmp(i)+tmp(i)-jm0(i) 
            sum1 = sum1+conjg(keep)*keep
            sum2 = sum2+conjg(keep)*jm1(i)
            jm0(i) = jm1(i)
            jm1(i) = keep
        end do
        moments(m+m) = moments(m+m)+real(sum2)
        moments(m+m+1) = moments(m+m+1)+real(sum1)
    end do
	!tmp = NULL
	!jm0 = NULL
	!jm1 = NULL
    deallocate(tmp,jm0,jm1)
end subroutine pem_diagonal_element
!********************************************************************

subroutine pem_sparse_product(rank, nonzeros, values, colind, rowptr, dest, src)
    implicit none
    integer :: rank, nonzeros, colind(nonzeros), rowptr(rank+1), i, j, k
    complex*16 :: values(nonzeros), dest(rank), src(rank)

    dest = (0.0,0.0)
    do j = 1, rank
        do k = rowptr(j), rowptr(j+1)-1
            i = colind(k)
            dest(j) = dest(j)+src(i)*values(k)
        end do
    end do
end subroutine pem_sparse_product

