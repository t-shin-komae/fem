module utils
    contains
        subroutine print_mat(mat)
            implicit none
            DOUBLE PRECISION , INTENT(IN) :: mat(:,:)
            INTEGER i,j
            do i = 1,ubound(mat,1)
                write (*,"(100f10.3) ") (mat(i,j),j=1,ubound(mat,1))
            end do
        end subroutine print_mat
        subroutine write_mat(mat,device)
            implicit none
            INTEGER , INTENT(IN) :: device
            DOUBLE PRECISION , INTENT(IN) :: mat(:,:)
            INTEGER i,j
            do i = 1,ubound(mat,1)
                write (device,"(100f10.3) ") (mat(i,j),j=1,ubound(mat,1))
            end do
        end subroutine write_mat
end module utils
