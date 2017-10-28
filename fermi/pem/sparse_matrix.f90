!Wirte by Liao Yuanda
!####################################################################

subroutine pem_sparse( H, values, colind, rowptr, rank, nonzeros)
	use globle_vab
	implicit none
	integer :: rank, nonzeros, colind(nonzeros), rowptr(rank+1), i, j, k
	complex*16 :: H(rank, rank), values(nonzeros)

		k = 1
		do i = 1, rank
			H(i,i) = H(i,i)-b
			rowptr(i) = k
			do j = 1, rank
				if ( H(i,j) /= (0.0,0.0)) then
					H(i,j) = H(i,j)/a
                    values(k) = H(i,j)
                    colind(k) = j
					k = k+1
                end if 
            end do
        end do

        rowptr(rank+1) = k
		write(*,*) k-1, nonzeros
end subroutine pem_sparse

subroutine diag_position(p, p_p, p_pN, pN_p, pN_pN, L)
    implicit none
    integer :: p, p_p, pN_p, p_pN, pN_pN, half_rank, L, nonzeros_perrow, i, j
    half_rank = L*L
    nonzeros_perrow = 8
	i = int((p-1)/L)
	j = mod((p-1),L)
    if (i == 0 .AND. j == 0) then
        p_p = 1
        p_PN = nonzeros_perrow
        pN_p = half_rank*nonzeros_perrow+1
        pN_PN = pN_p+1
    else if(i == 0 .AND. j /= 0 .AND. j /= L-1) then
        p_p = (p-1)*nonzeros_perrow+2
        p_PN = p*nonzeros_perrow
        pN_p = (p-1+half_rank)*nonzeros_perrow+1
        pN_pN = pN_p+2
    else if( (i == 0 .AND. j == L-1) .OR. (i /= 0 .AND. i /= L-1 .AND. j == 0) ) then
        p_p = (p-1)*nonzeros_perrow+3
        p_PN = p*nonzeros_perrow
        pN_p = (p-1+half_rank)*nonzeros_perrow+1
        pN_pN = pN_p+3
    else if(i /= 0 .AND. i /= L-1 .AND. j /= 0 .AND. j /= L-1) then
        p_p = (p-1)*nonzeros_perrow+4
        p_PN = p*nonzeros_perrow
        pN_p = (p-1+half_rank)*nonzeros_perrow+1
        pN_pN = pN_p+4
    else if( (i /= 0 .AND. i /= L-1 .AND. j == L-1) .OR. (i == L-1 .AND. j == 0) ) then
        p_p = (p-1)*nonzeros_perrow+5
        p_PN = p*nonzeros_perrow
        pN_p = (p-1+half_rank)*nonzeros_perrow+1
        pN_pN = pN_p+5
    else if( i == L-1 .AND. j /= 0 .AND. j /= L-1 ) then
        p_p = (p-1)*nonzeros_perrow+6
        p_PN = p*nonzeros_perrow
        pN_p = (p-1+half_rank)*nonzeros_perrow+1
        pN_pN = pN_p+6
    else if( i == L-1 .AND. j == L-1 ) then
        p_p = (p-1)*nonzeros_perrow+7
        p_PN = p*nonzeros_perrow
        pN_p = (p-1+half_rank)*nonzeros_perrow+1
        pN_pN = pN_p+7
    end if
end subroutine diag_position
