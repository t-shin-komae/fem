module mesh
    contains
        subroutine dirichlet(A,b,nth,alpha)
            implicit none
            DOUBLE PRECISION , INTENT(INOUT) :: A(:,:),b(:)
            DOUBLE PRECISION , INTENT(IN) :: alpha
            INTEGER , INTENT(IN) :: nth
            b(:) = b(:) - A(:,nth) * alpha
            b(nth) = alpha
            A(:,nth) = 0d0
            A(nth,:) = 0d0
            A(nth,nth) = 1.0d0
        end subroutine dirichlet
end module