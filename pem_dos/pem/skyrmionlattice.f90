!Write by Liao Yuanda
!
!###############################################################

    subroutine calcu_H( L, N, M, H, t )

        implicit none
        integer :: L, N, M, i, j
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl
        complex*16  :: H(M, M)
		real*8 :: t

        H = 0.0

        do i = 0, L-1, 1
            do j = 0, L-1, 1
                p = i*L+j+1
                call triangular_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )

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

program skyrmionlattice
  	use globle_vab
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !##########################################
    integer, parameter :: L = 60, N = L*L, M = 2*N, nonzeros = 7*M, cutoff = 100, num_x = 100
    real*8, parameter :: t = 1.0
    integer, dimension(nonzeros) :: colind
    integer, dimension(M+1) :: rowptr
    complex*16, dimension(nonzeros) :: values
    complex*16, allocatable :: H(:,:)
	real*8, dimension(cutoff) :: miun_pem
	!real*8, dimension(num_x) :: stateden_pem_integral
    !real*8, dimension(cutoff) :: miun_tpem, miun_pem, miun_ev
    !real*8, dimension(M) :: eigv
    !real*8, dimension(num_x) :: stateden_ev, stateden_pem
    integer :: i, j, interval
    real*8 :: weight, t1, t2

    Emax = 3*t
    Emin = -6*t
    a = (Emax-Emin)/2.0
    b = (Emax+Emin)/2.0
    !##########################
    call cpu_time(t1)
    write(*,*) 'Begin at:', t1
    !######################
	!interval = 100

	!open(unit = 15, file = 'eigv.in')
    open(unit = 14, file = 'miun_pem.out')
	!open(unit = 13, file = 'dos_ev.out')
	!open(unit = 17, file = 'dos_pem.out')
	!open(unit = 18, file = 'dos_pem_integral.out')
	!do i = 1, M
	!	read(15,*) eigv(i)
	!end do
	!close(15)
    allocate(H(M,M))
    call calcu_H( L, N, M, H, t )
    call pem_sparse( H, values, colind, rowptr, M, nonzeros)
    deallocate(H)
	!write(*,*) 1
	!call pem_statedensity( num_x, cutoff, M, nonzeros, val1, val2, colind, rowptr, stateden_pem )
	!call pem_statedensity_integral( num_x, interval, cutoff, M, nonzeros, val1, val2, colind, rowptr, stateden_pem_integral)
	!call pem_particle_density( cutoff, M, nonzeros, val1, val2, colind, rowptr, beta, miu, a, b, weight )
	call pem_calculate_moments(cutoff, M, nonzeros, values, colind, rowptr, miun_pem)
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
