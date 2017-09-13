!Write by Liao Yuanda
!
!###############################################################

module dou_ex
contains

    subroutine cau_Hf( L, N, M, Hf, t )
        !###############################################
        !(i,j) to (i, j_new); (i,j) to (i_new,j); N = i+(j-1)*L
        !###############################################
        implicit none
        integer :: L, N, M, i, j, k, i_minus, i_plus, j_minus, j_plus, k_minus, k_plus, N_ij, N_i_plus, N_i_minus, N_j_minus, N_j_plus, N_k_minus, N_k_plus
        real*8 :: t, Hf(M, M)
        do i = 1, M
            do j = 1, M
                Hf(i,j) = 0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                do k = 0, L-1, 1
                !period of boundary 
                i_plus = mod(i+1, L)
                j_plus = mod(j+1, L)
                k_plus = mod(k+1, L)
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
                if (k==0) then
                    k_minus = L-1
                else
                    k_minus = k-1
                end if
                N_ij = i*L*L+j*L+k+1
                N_i_minus = i_minus*L*L+j*L+k+1
                N_i_plus = i_plus*L*L+j*L+k+1
                N_j_minus = i*L*L+j_minus*L+k+1
                N_j_plus = i*L*L+j_plus*L+k+1
                N_k_minus = i*L*L+j*L+k_minus+1
                N_k_plus = i*L*L+j*L+k_plus+1

                Hf(N_ij, N_i_minus) = -t
                Hf(N_ij, N_i_plus) = -t
                Hf(N_ij, N_j_minus) = -t
                Hf(N_ij, N_j_plus) = -t

                Hf(N_ij+N, N_i_minus+N) = -t
                Hf(N_ij+N, N_i_plus+N) = -t
                Hf(N_ij+N, N_j_minus+N) = -t
                Hf(N_ij+N, N_j_plus+N) = -t
                if ( k == 0 ) then
                    Hf(N_ij, N_k_minus) = t
                    Hf(N_ij+N, N_k_minus+N) = t
                elseif ( k == L-1 ) then
                    Hf(N_ij, N_k_plus) = t
                    Hf(N_ij+N, N_k_plus+N) = t
                else
                    Hf(N_ij, N_k_minus) = -t
                    Hf(N_ij, N_k_plus) = -t
                    Hf(N_ij+N, N_k_minus+N) = -t
                    Hf(N_ij+N, N_k_plus+N) = -t
                end if
               
                end do
            end do
        end do
    end subroutine cau_Hf

    subroutine cau_H( N, M, Hf, H, J_para, fai, theta )
        implicit none
        integer :: N, M, i, j, nonzeros, zeros
        real*8 :: J_para, fai(N), theta(N), Hf(M,M)
        complex*16 :: H(M, M)
        nonzeros = 0
        zeros = 0
        do i = 1, M
            do j = 1, M
                H(i,j) = (0.0, 0.0)
            end do
        end do

        do i = 1, N 
            H(i, i) = -J_para*cos(theta(i)) 
            H(i+N, i+N) = J_para*cos(theta(i)) 
            H(i+N, i) = -J_para*sin(theta(i))*exp( cmplx(0.0, fai(i)) ) 
            H(i, i+N) = -J_para*sin(theta(i))*exp( cmplx(0.0, fai(i)*(-1)) ) 
        end do
        do i = 1, M
            do j = 1, M
                H(i,j) = Hf(i,j)+H(i,j)
            end do
        end do

    end subroutine cau_H
    
    subroutine cau_temp_sparse( H, val1, val2, colind, rowptr , M, a, b, nonzeros)
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

    end subroutine
    
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


program double_exchange
    use dou_ex 
    implicit none
    external ctofortran
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_para, the parameter if H
    !##########################################
    integer, parameter :: L = 4, N = 64, M = 128, num_MC = 10**5, nonzeros = 8*M
    real*8, parameter :: t = 1.0, PI = 3.141592654
    integer :: i, j, k, ii, p, num_aver, d, di, dd, ij
    real*8 :: sx, sy, sz, ran, J_para, beta, miu, Tc, x, ln_Weight, ln_Weight_new, ctofortran, Emax, Emin, a, b
    integer, dimension(nonzeros) :: colind
    integer, dimension(M+1) :: rowptr
    real*8, dimension(2,N) :: fai, theta
    real*8, dimension(nonzeros) :: val1, val2
    real*8, dimension(int(num_MC/N)) :: ss, s_M
    real*8, dimension(M,M) :: Hf
    complex*16, dimension(M,M) :: H
    J_para = 8*t
    miu = -8*t
    d = 1
    call initran(1)
    call cau_Hf( L, N, M, Hf, t )
    open(unit=13, file='para.in', status='old')
    open(unit=12, file='out.dat')
    read(13,*) Tc 
    Tc = Tc*t*0.03
    beta = 1/Tc
    Emax = 2*L*t+J_para
    Emin = -2*L*t-J_para
    a = (Emax-Emin)/2.0
    b = (Emax+Emin)/2.0
    ij = 0
    do i = 1, N
        x = ran()
        fai(d,i) = x*2*PI
        x = ran()
        theta(d,i) = x*PI 
    end do
    call cau_H( N, M, Hf, H, J_para, fai(d,:), theta(d,:))
    !write(*,*) nonzeros
    call cau_temp_sparse( H, val1, val2, colind, rowptr , M, a, b, nonzeros)
    ln_Weight =  ctofortran(M, nonzeros, val1, val2, colind, rowptr, beta, miu, a, b)


    do i = 1, num_MC 
        dd = di(d)
        x = ran()
        p = int(x*N)+1 
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
        call cau_H( N, M, Hf, H, J_para, fai(dd,:), theta(dd,:))
        call cau_temp_sparse( H, val1, val2, colind, rowptr , M, a, b, nonzeros)
        ln_Weight_new = ctofortran(M, nonzeros, val1, val2, colind, rowptr, beta, miu, a, b)

        x = ran()
        !write(*,*) x, '<?', exp( ln_Weight_new-ln_Weight )
        !write(*,*) ln_Weight, ln_Weight_new
        if ( x > exp( ln_Weight_new-ln_Weight ) ) then
            fai(dd,p) = fai(d,p)
            theta(dd, p) = theta(d, p)
        else
            ln_Weight = ln_Weight_new
        end if
    if ( mod(i,N)==0 ) then
        ij = ij+1
        ss(ij) = 0
        sx = 0
        sy = 0 
        sz = 0
        s_M(ij) = 0
        do j = 1, N
            sx = sx+sin(theta(d,j))*cos(fai(d,j))
            sy = sy+sin(theta(d,j))*sin(fai(d,j))
            sz = sz+cos(theta(d,j))
            do k = 1, N
                ss(ij) = ss(ij)+sin(theta(d,j))*cos(fai(d,j))*sin(theta(d,k))*cos(fai(d,k))
                ss(ij) = ss(ij)+sin(theta(d,j))*sin(fai(d,j))*sin(theta(d,k))*sin(fai(d,k))
                ss(ij) = ss(ij)+cos(theta(d,j))*cos(theta(d,k))
            end do
        end do
        s_M(ij) = sqrt(sx**2+sy**2+sz**2)/N
        ss(ij) = ss(ij)/N/N
    end if
        d = dd
        !write(*,*) i/(num_MC/100)
    end do

    sx = 0
    sy = 0
    num_aver = 1000
    do i = 1, num_aver 
        sx = sx+ss(ij-num_aver+i)
        sy = sy+s_M(ij-num_aver+i)
    end do
    sx = sx/num_aver
    sy = sy/num_aver
    write(12,*) Tc, sx, sy

    write(*,*) 'end'
end program double_exchange
