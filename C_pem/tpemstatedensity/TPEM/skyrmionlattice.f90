!Write by Liao Yuanda
!
!###############################################################

module dou_ex

contains

    subroutine calcu_p_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )
        implicit none
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L
        integer :: i, j, i_pl, j_pl, i_mi, j_mi

        i = int((p-1)/L)
        j = mod(p-1,L)
        i_pl = mod(i+1,L)
        j_pl = mod(j+1,L)
        if (i==0) then
            i_mi = L-1
        else
            i_mi = i-1
        end if
        if (j==0) then
            j_mi = L-1
        else
            j_mi = j-1
        end if

        p_i_pl = i_pl*L+j+1
        p_i_mi = i_mi*L+j+1
        p_j_pl = i*L+j_pl+1
        p_j_mi = i*L+j_mi+1
        p_i_pl_j_mi = i_pl*L+j_mi+1
        p_i_mi_j_pl = i_mi*L+j_pl+1

    end subroutine calcu_p_nearest 

    subroutine calcu_H( L, N, M, H, t )

        implicit none
        integer :: L, N, M, i, j
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl
        complex*16  :: H(M, M)
		real*8 :: t

        do i = 1, M
            do j = 1, M
                H(i,j) = 0.0
            end do
        end do

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                p = i*L+j+1
                call calcu_p_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )

                H(p, p_i_pl) = -t
                H(p, p_i_mi) = -t
                H(p, p_j_pl) = -t
                H(p, p_j_mi) = -t
                H(p, p_i_pl_j_mi) = -t
                H(p, p_i_mi_j_pl) = -t

                H(p+N, p_i_pl+N) = -t
                H(p+N, p_i_mi+N) = -t
                H(p+N, p_j_pl+N) = -t
                H(p+N, p_j_mi+N) = -t
                H(p+N, p_i_pl_j_mi+N) = -t
                H(p+N, p_i_mi_j_pl+N) = -t
            end do
        end do
    end subroutine calcu_H

    subroutine pem_sparse( H, val1, val2, colind, rowptr , M, nonzeros)
            implicit none
            integer :: M, nonzeros, colind(nonzeros), rowptr(M+1), i, j, k, num_row
            complex*16 :: H(M, M)
            real*8 :: val1(nonzeros), val2(nonzeros)
            k = 0
            do i = 1, M
                rowptr(i) = k

                do j = 1, M
                    if ( H(i,j) /= (0.0,0.0)) then
                        k = k+1

                        val1(k) = real(H(i,j))
                        val2(k) = aimag(H(i,j))
                        colind(k) = j-1
                    end if 
                end do
            end do
            rowptr(M+1) = nonzeros
			write(*,*) k, nonzeros
    end subroutine pem_sparse
	
	subroutine pem_rescale( H, a, b, M)
		implicit none
		integer :: M, i, j
		complex*16 :: H(M,M)
		real*8 :: a, b
		do i = 1,M
			H(i,i) = H(i,i)-b
        	do j = 1,M
				H(i,j) = H(i,j)/a
        	end do
    	end do
	end subroutine pem_rescale
end module dou_ex

program skyrmionlattice
    use dou_ex
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_para, the parameter if H
    !##########################################
    integer, parameter :: L = 60, N = L*L, M = 2*N, nonzeros = 7*M, cutoff = 800, num_x = 100
    real*8, parameter :: t = 1.0
    integer, dimension(nonzeros) :: colind
    integer, dimension(M+1) :: rowptr
    real*8, dimension(nonzeros) :: val1, val2
	real*8, dimension(cutoff) :: miun_pem
	!real*8, dimension(num_x) :: stateden_pem_integral
    !real*8, dimension(cutoff) :: miun_tpem, miun_pem, miun_ev
    !real*8, dimension(M) :: eigv
    !real*8, dimension(num_x) :: stateden_ev, stateden_pem
    integer :: i, j, interval
    real*8 :: Emax, Emin, a, b, weight, beta, miu, t1, t2
    complex*16, dimension(M,M) :: H
    Emax = 3*t
    Emin = -6*t
    a = (Emax-Emin)/2.0
    b = (Emax+Emin)/2.0
	beta = 10.0*t
	miu = 0.7424*t
	!interval = 100
    !##########################
    call cpu_time(t1)
    write(*,*) 'Begin at:', t1
    !######################
	!open(unit = 15, file = 'eigv.in')
    open(unit = 14, file = 'miun_pem.out')
	!open(unit = 13, file = 'dos_ev.out')
	!open(unit = 17, file = 'dos_pem.out')
	!open(unit = 18, file = 'dos_pem_integral.out')
	!do i = 1, M
	!	read(15,*) eigv(i)
	!end do
	!close(15)
    call calcu_H( L, N, M, H, t )
	call pem_rescale( H, a, b, M)
    call pem_sparse( H, val1, val2, colind, rowptr , M, nonzeros)

	!call pem_statedensity( num_x, cutoff, M, nonzeros, val1, val2, colind, rowptr, stateden_pem )
	!call pem_statedensity_integral( num_x, interval, cutoff, M, nonzeros, val1, val2, colind, rowptr, stateden_pem_integral)
	!call pem_particle_density( cutoff, M, nonzeros, val1, val2, colind, rowptr, beta, miu, a, b, weight )
	call pem_moments( cutoff, M, nonzeros, val1, val2, colind, rowptr, miun_pem)
	!call pem_logweight( cutoff, M, nonzeros, val1, val2, colind, rowptr, 0.0, 0.0, a, b, weight )
	!call pem_test( cutoff, M, nonzeros, val1, val2, colind, rowptr, miun_tpem, miun_pem, miun_ev, eigv)
	!call pem_ev_dos( cutoff, num_x, M, eigv, stateden_ev)
	!write(*,*) weight
    !do i = 1, num_x
		!write(13,*) stateden_ev(i)
		!write(17,*) stateden_pem(i)
		!write(18,*) stateden_pem_integral(i)
    !end do
	!write(*,*) weight
    do i = 1, cutoff
		!write(14,*) miun_tpem(i),miun_pem(i),miun_ev(i)
		write(14,*) miun_pem(i)
    end do
    !close(14)
    !close(13)
    !close(17)
    !close(18)
    !##########################
    call cpu_time(t2)
    write(*,*) 'End at:', t2
    !######################
    write(*,*) 'Timecost:', t2-t1
end program skyrmionlattice
