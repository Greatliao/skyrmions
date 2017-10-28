
	!Triangular lattice 2D
    subroutine triangular_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_i_pl_j_mi, p_i_mi_j_pl, L )
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

    end subroutine Triangular_nearest

	!cubic lattice 
    subroutine cubic_nearest( p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_k_pl, p_k_mi, L )
        implicit none
        integer :: p, p_i_pl, p_i_mi, p_j_pl, p_j_mi, p_k_pl, p_k_mi, L
        integer :: i, j, k, i_pl, j_pl, k_pl, i_mi, j_mi, k_mi

        i = (p-1)/L/L
        k = mod(p-1,L)
		j = (p-1-i*L*L)/L
        i_pl = mod(i+1,L)
        j_pl = mod(j+1,L)
		k_pl = mod(k+1,L)
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
        if (k==0) then
            k_mi = L-1
        else
            k_mi = k-1
        end if

        p_i_pl = i_pl*L*L+j*L+k+1
        p_i_mi = i_mi*L*L+j*L+k+1
        p_j_pl = i*L*L+j_pl*L+k+1
        p_j_mi = i*L*L+j_mi*L+k+1
        p_k_pl = i*L*L+j*L+k_pl+1
        p_k_mi = i*L*L+j*L+k_mi+1

    end subroutine cubic_nearest
