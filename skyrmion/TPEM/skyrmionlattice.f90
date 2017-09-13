!Write by Liao Yuanda
!
!###############################################################

module dou_ex
contains

    subroutine calcu_Hf( L, N, M, Hf, t )
        !###############################################
        !(i,j) to (i, j_new); (i,j) to (i_new,j); N = i+(j-1)*L
        !###############################################
        implicit none
        integer :: L, N, M, i, j, i_minus, i_plus, j_minus, j_plus, N_ij, N_i_plus, N_i_minus, N_j_minus, N_j_plus, N_i_minus_j_plus, N_i_plus_j_minus
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) = 0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                !period of boundary 
                i_plus = mod(i+1, L)
                j_plus = mod(j+1, L)
                if (i==0) then
                    i_minus = L-1
                else
                    i_minus = i-1
                end if 
                if (j==0) then
                    j_minus = L-1
                else
                    j_minus = j-1
                end if 
                N_ij = i*L+j+1
                N_i_minus = i_minus*L+j+1
                N_i_plus = i_plus*L+j+1
                N_j_minus = i*L+j_minus+1
                N_j_plus = i*L+j_plus+1
                N_i_minus_j_plus = i_minus*L+j_plus+1
                N_i_plus_j_minus = i_plus*L+j_minus+1

                Hf(N_ij, N_i_minus) = -t
                Hf(N_ij, N_i_plus) = -t
                Hf(N_ij, N_j_minus) = -t
                Hf(N_ij, N_j_plus) = -t
                Hf(N_ij, N_i_plus_j_minus) = -t
                Hf(N_ij, N_i_minus_j_plus) = -t

                Hf(N_ij+N, N_i_minus+N) = -t
                Hf(N_ij+N, N_i_plus+N) = -t
                Hf(N_ij+N, N_j_minus+N) = -t
                Hf(N_ij+N, N_j_plus+N) = -t
                Hf(N_ij+N, N_i_plus_j_minus+N) = -t
                Hf(N_ij+N, N_i_minus_j_plus+N) = -t
               
            end do
        end do
    end subroutine calcu_Hf

    subroutine calcu_H( N, M, Hf, H, J_k, fai, theta )
        implicit none
        integer :: N, M, i, j
        real*8 :: J_k, fai(N), theta(N), Hf(M,M)
        complex*16 :: H(M, M)
        do i = 1, M
            do j = 1, M
                H(i,j) = (0.0, 0.0)
            end do
        end do

        do i = 1, N 
            H(i, i) = -J_k*cos(theta(i)) 
            H(i+N, i+N) = J_k*cos(theta(i)) 
            H(i+N, i) = -J_k*sin(theta(i))*exp( cmplx(0.0, fai(i)) ) 
            H(i, i+N) = -J_k*sin(theta(i))*exp( cmplx(0.0, fai(i)*(-1)) ) 
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = Hf(i,j)+H(i,j)
            end do
        end do

    end subroutine calcu_H

    subroutine calcu_V( L, N, J_H, lambda, fai, theta, V )
        implicit none
        integer :: L, N, i, j, i_minus, i_plus, j_minus, j_plus, N_i_minus, N_i_plus, N_j_minus, N_j_plus, N_ij, N_i_plus_j_minus, N_i_minus_j_plus
        real*8 :: J_H, lambda, fai(N), theta(N), V
        V = 0
        do i = 0, L-1, 1
            do j = 0, L-1, 1
                !period of boundary 
                i_plus = mod(i+1, L)
                j_plus = mod(j+1, L)
                if (i==0) then
                    i_minus = L-1
                else
                    i_minus = i-1
                end if 
                if (j==0) then
                    j_minus = L-1
                else
                    j_minus = j-1
                end if 
                N_ij = i*L+j+1
                N_i_minus = i_minus*L+j+1
                N_i_plus = i_plus*L+j+1
                N_j_minus = i*L+j_minus+1
                N_j_plus = i*L+j_plus+1
                N_i_minus_j_plus = i_minus*L+j_plus+1
                N_i_plus_j_minus = i_plus*L+j_minus+1

                V = V-1.0/2*J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_i_minus))*cos(fai(N_i_minus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_i_plus))*cos(fai(N_i_plus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_j_minus))*cos(fai(N_j_minus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_j_plus))*cos(fai(N_j_plus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_i_plus_j_minus))*cos(fai(N_i_plus_j_minus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*cos(fai(N_ij))*sin(theta(N_i_minus_j_plus))*cos(fai(N_i_minus_j_plus))

                V = V-1.0/2*J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_i_minus))*sin(fai(N_i_minus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_i_plus))*sin(fai(N_i_plus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_j_minus))*sin(fai(N_j_minus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_j_plus))*sin(fai(N_j_plus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_i_plus_j_minus))*sin(fai(N_i_plus_j_minus))
                V = V-1.0/2*J_H*sin(theta(N_ij))*sin(fai(N_ij))*sin(theta(N_i_minus_j_plus))*sin(fai(N_i_minus_j_plus))

                V = V-1.0/2*J_H*cos(theta(N_ij))*cos(theta(N_i_minus))
                V = V-1.0/2*J_H*cos(theta(N_ij))*cos(theta(N_i_plus))
                V = V-1.0/2*J_H*cos(theta(N_ij))*cos(theta(N_j_minus))
                V = V-1.0/2*J_H*cos(theta(N_ij))*cos(theta(N_j_plus))
                V = V-1.0/2*J_H*cos(theta(N_ij))*cos(theta(N_i_plus_j_minus))
                V = V-1.0/2*J_H*cos(theta(N_ij))*cos(theta(N_i_minus_j_plus))

                V = V+1.0/2*lambda*(-1.0)*( sin(theta(N_i_minus))*sin(fai(N_i_minus))*cos(theta(N_ij))-cos(theta(N_i_minus))&
                    *sin(theta(N_ij))*sin(fai(N_ij)) )
                V = V+1.0/2*lambda*( sin(theta(N_i_plus))*sin(fai(N_i_plus))*cos(theta(N_ij))-cos(theta(N_i_plus))&
                    *sin(theta(N_ij))*sin(fai(N_ij)) )
                V = V+1.0/2*lambda*(-1.0)*( (1.0/2)*( sin(theta(N_j_minus))*sin(fai(N_j_minus))*cos(theta(N_ij))&
                    -cos(theta(N_j_minus))*sin(theta(N_ij))*sin(fai(N_ij)) ) + sqrt(3.0)/2*( cos(theta(N_j_minus))&
                    *sin(theta(N_ij))*cos(fai(N_ij))-sin(theta(N_j_minus))*cos(fai(N_j_minus))*cos(theta(N_ij)) ) )
                V = V+1.0/2*lambda*( (1.0/2)*( sin(theta(N_j_plus))*sin(fai(N_j_plus))*cos(theta(N_ij))&
                    -cos(theta(N_j_plus))*sin(theta(N_ij))*sin(fai(N_ij)) ) &
                    + sqrt(3.0)/2*( cos(theta(N_j_plus))*sin(theta(N_ij))*cos(fai(N_ij))-sin(theta(N_j_plus))*cos(fai(N_j_plus))*cos(theta(N_ij)) ) )
                V = V+1.0/2*lambda*( sin(theta(N_i_plus_j_minus))*sin(fai(N_i_plus_j_minus))*cos(theta(N_ij))&
                    -cos(theta(N_i_plus_j_minus))*sin(theta(N_ij))*sin(fai(N_ij)) ) &
                    +1.0/2*lambda*(-1.0)*( (1.0/2)*( sin(theta(N_i_plus_j_minus))*sin(fai(N_i_plus_j_minus))*cos(theta(N_ij))-cos(theta(N_i_plus_j_minus))*sin(theta(N_ij))*sin(fai(N_ij)) ) &
                    + sqrt(3.0)/2*( cos(theta(N_i_plus_j_minus))*sin(theta(N_ij))*cos(fai(N_ij))-sin(theta(N_i_plus_j_minus))*cos(fai(N_i_plus_j_minus))*cos(theta(N_ij)) ) )
                V = V+1.0/2*lambda*(-1.0)*( sin(theta(N_i_minus_j_plus))*sin(fai(N_i_minus_j_plus))*cos(theta(N_ij))&
                    -cos(theta(N_i_minus_j_plus))*sin(theta(N_ij))*sin(fai(N_ij)) ) &
                    +1.0/2*lambda*( (1.0/2)*( sin(theta(N_i_minus_j_plus))*sin(fai(N_i_minus_j_plus))*cos(theta(N_ij))-cos(theta(N_i_minus_j_plus))&
                    *sin(theta(N_ij))*sin(fai(N_ij)) ) + sqrt(3.0)/2*( cos(theta(N_i_minus_j_plus))*sin(theta&
                    (N_ij))*cos(fai(N_ij))-sin(theta(N_i_minus_j_plus))*cos(fai(N_i_minus_j_plus))*cos(theta(N_ij)) ) )

            end do
        end do


    end subroutine calcu_V

    subroutine calcu_deltaV(L, N, B_z, theta, deltaV)
        implicit none
        integer :: L, i, j, N_ij, N
        real*8 :: B_z, theta(N), deltaV
        deltaV = 0
        do i = 0, L-1, 1
            do j = 0, L-1, 1
                N_ij = i*L+j+1
                deltaV = deltaV+B_z*cos(theta(N_ij))
            end do
        end do

    end subroutine calcu_deltaV
    
    subroutine calcu_temp_sparse( H, val1, val2, colind, rowptr , M, a, b, nonzeros)
            implicit none
            integer :: M, nonzeros, colind(nonzeros), rowptr(M+1), i, j, k, num_row
            complex*16 :: H(M, M)
            real*8 :: val1(nonzeros), val2(nonzeros), a, b
            k = 0
            do i = 1, M
                rowptr(i) = k
                num_row = 0
                do j = 1, M
                    if ( H(i,j) /= (0.0,0.0)) then
                        k = k+1
                        num_row = num_row+1
                        val1(k) = (real(H(i,j))-b)/a
                        val2(k) = aimag(H(i,j))/a
                        colind(k) = j-1
                    end if 
                end do
            end do
            rowptr(M+1) = nonzeros

    end subroutine calcu_temp_sparse
    
end module dou_ex

real(8) function ran()
!----------------------------------------------!
! 64-bit congruental generator                 !
! iran64=oran64*2862933555777941757+1013904243 !
!----------------------------------------------!
    implicit none

    real(8)    :: dmu64
    integer(8) :: ran64,mul64,add64
    common/bran64/dmu64,ran64,mul64,add64

    ran64=ran64*mul64+add64
    ran=0.5d0+dmu64*dble(ran64)

end function ran
!----------------!

!---------------------!
subroutine initran(w)
!---------------------!
    implicit none

    integer(8) :: irmax
    integer(4) :: w,nb,b

    real(8)    :: dmu64
    integer(8) :: ran64,mul64,add64
    common/bran64/dmu64,ran64,mul64,add64
                        
    irmax=2_8**31
    irmax=2*(irmax**2-1)+1
    mul64=2862933555777941757_8
    add64=1013904243
    dmu64=0.5d0/dble(irmax)

    open(10,file='seed.in',status='old')

    read(10,*)ran64
    close(10)
    if (w.ne.0) then
        open(10,file='seed.in',status='unknown')
        write(10,*) abs((ran64*mul64)/5+5265361)
        close(10)
    end if

end subroutine initran
!----------------------!

!function to less dimension
integer function di(d)
    integer :: d
    if (d == 1) then
        di = 2
    end if
    if ( d == 2) then
        di =1
    end if
end function


program skyrmionlattice
    use dou_ex 
    implicit none
    external ctofortran
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_k, the parameter if H
    !##########################################
    integer, parameter :: L = 21, N = L*L, M = 2*N, num_T_para = 20, num_MC = 10**6, nonzeros = 8*M
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, ii, p, d, di, dd, ij
    real*8 :: ran, J_k, beta, miu, Tc, x, ln_Weight, ln_Weight_new, ln_Weight_fermi, ctofortran, Emax, Emin, a, b
    real*8 :: V, deltaV, ln_Weight_fermi_new, J_H, lambda, B_z, pa
    integer, dimension(nonzeros) :: colind
    integer, dimension(M+1) :: rowptr
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(nonzeros) :: val1, val2
    real*8, dimension(M,M) :: Hf
    real*8, dimension(num_T_para) :: T_para
    complex*16, dimension(M,M) :: H
    character( len = 4 ) :: cha
    !open( unit =12, file = 'para.in')
    !read(12,*) pa
    J_k = 8*t
    miu = -8*t
    d = 1
    J_H = 1.0*t
    lambda = 13/21*J_H*2*PI
    !B_z = lambda**2/J_H
    B_z = 0

    T_para = (/1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2 ,0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, &
    0.009/)
    call initran(1)
    call calcu_Hf( L, N, M, Hf, t )

    Emax = 6*t+J_k
    Emin = -6*t-J_k
    a = (Emax-Emin)/2.0
    b = (Emax+Emin)/2.0
    do i = 1, N
        x = ran()
        fai(d,i) = x*2*PI
        x = ran()
        theta(d,i) = x*PI
    end do
do ii = 1, num_T_para
    !ij = 0
    write(cha,'(i4)') ii
    open( unit = ii+20, file = 'lattice'//adjustL( trim(cha) )//'.out' )
    Tc = T_para(ii)*t
    beta = 1/Tc
    call calcu_H( N, M, Hf, H, J_k, fai(d,:), theta(d,:))
    !write(*,*) nonzeros
    call calcu_temp_sparse( H, val1, val2, colind, rowptr , M, a, b, nonzeros)
    call calcu_V( L, N, J_H, lambda, fai(d,:), theta(d,:), V )
    call calcu_deltaV(L, N, B_z, theta(d,:), deltaV)
    ln_Weight_fermi =  ctofortran(M, nonzeros, val1, val2, colind, rowptr, beta, miu, a, b)
    ln_Weight = ln_Weight_fermi-beta*(V+deltaV)
    !p = 0
    do i = 1, num_MC 
        dd = di(d)
        x = ran()
        p = int(x*N)+1
        !p =  p+1
        !Write(*,*) p
        fai(dd,:) = fai(d,:)
        theta(dd,:) = theta(d,:)
        x = ran()
        fai(dd, p) = fai(d,p)+(x-0.5)*2*2*PI/12
        if ( fai(dd,p) < 0 ) then
            fai(dd,p) = fai(dd,p)+2*PI
        end if
        if ( fai(dd,p) > 2*PI ) then
            fai(dd,p) = fai(dd,p)-2*PI
        end if

        x = ran()
        theta(dd, p) = theta(d,p)+(x-0.5)*2*PI/6
        if ( theta(dd,p) < 0 ) then
            theta(dd,p) = theta(dd,p)+PI
        end if
        if (theta(dd,p) > PI ) then
            theta(dd,p) = theta(dd,p)-PI
        end if
        call calcu_H( N, M, Hf, H, J_k, fai(dd,:), theta(dd,:))
        call calcu_temp_sparse( H, val1, val2, colind, rowptr , M, a, b, nonzeros)
        call calcu_V( L, N, J_H, lambda, fai(dd,:), theta(dd,:), V )
        call calcu_deltaV(L, N, B_z, theta(dd,:), deltaV)
        ln_Weight_fermi_new = ctofortran(M, nonzeros, val1, val2, colind, rowptr, beta, miu, a, b)
        ln_Weight_new = ln_Weight_fermi_new-beta*(V+deltaV)

        x = ran()
        !write(*,*) x, '<?', exp( ln_Weight_new-ln_Weight )
        !write(*,*) ln_Weight, ln_Weight_new
        if ( x > exp( ln_Weight_new-ln_Weight ) ) then
            fai(dd,p) = fai(d,p)
            theta(dd, p) = theta(d, p)
        else
            ln_Weight = ln_Weight_new
        end if
        d = dd
        !write(*,*) i/(num_MC/100)
        if ( i > num_MC-10*N ) then

            if ( mod(i+10*N-num_MC, 100) == 0 ) then
                !ij = ij+1
                !write(ii+20,*) ij
                do j = 1, N
                    write(ii+20,*) j, theta(d,j), fai(d,j)
                end do
            end if
        end if
        !if ( p == N ) then
        !    p = 0
        !end if

    end do
    close(ii+20)
end do


    write(*,*) 'end'
end program skyrmionlattice
