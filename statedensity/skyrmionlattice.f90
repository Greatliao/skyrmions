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
        real*8 :: t, H(M, M)

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
                H(p, p_i_mi) = -t
                H(p, p_i_pl_j_mi) = -t
                H(p, p_i_mi_j_pl) = -t

                H(p+N, p_i_pl+N) = -t
                H(p+N, p_i_mi+N) = -t
                H(p+N, p_j_pl+N) = -t
                H(p+N, p_i_mi+N) = -t
                H(p+N, p_i_pl_j_mi+N) = -t
                H(p+N, p_i_mi_j_pl+N) = -t
            end do
        end do
    end subroutine calcu_H

end module dou_ex

program skyrmionlattice
    use dou_ex
    implicit none
    !#########################################
    !L, size of lattice
    !N = L*L, num of lattice
    !t, J_para, the parameter if H
    !##########################################
    integer, parameter :: L = 60, N = L*L, M = 2*N, LDZ = 1
    real*8, parameter :: t = 1.0
    integer :: i, j, info
    real*8, dimension(M) :: Eig
    real*8, dimension(3*M-2) :: Rwork
    real*8, dimension(M,M) :: H
    real*8, dimension(3*M) :: Work
    real*8, dimension(M*(M+1)/2) :: AP
    real*8, dimension(LDZ,M) :: Z

    call calcu_H( L, N, M, H, t )
    do i = 1,M
        do j = i,M
            AP(i+(j-1)*j/2) = H(i,j)
        end do
    end do
    call dspev( 'N', 'U', M, AP, Eig, Z, LDZ, work, info )
    if (info /= 0) then
    write(*,*) 'info_erro =', info
    end if
    Eig = Eig/t
    open(unit = 14, file = 'eig.out')
    do i = 1, M
        write(14,*) Eig(i)
    end do
    close(14)

    write(*,*) 'end'
end program skyrmionlattice
