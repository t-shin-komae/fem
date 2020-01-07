module two_dimentional
implicit none
type node
    INTEGER id
    DOUBLE PRECISION point(2)
    ! DOUBLE PRECISION f
end type
type triangle
    PRIVATE
    TYPE(node) :: nodes(6)
end type

contains
subroutine create_quad_matrix(K,f,n,triangles)
    implicit none
    INTEGER , INTENT(IN) :: n ! size of matrix A(n*n) and vector b
    DOUBLE PRECISION , INTENT(OUT) :: K(n,n),f(n)
    TYPE(triangle) :: triangles(:),triangle_i
    DOUBLE PRECISION :: points_i(6,2),Ke(6,6)
    INTEGER :: ids_i(6)
    INTEGER i,j,l

    K(:,:) = 0
    f(:) = 0

    do i = 1,ubound(triangles,1)
        triangle_i = triangles(i)
        ids_i = triangle_i%nodes%id
        do l = 1,6
            points_i(l,:) = triangle_i%nodes(l)%point
        end do
        Ke = create_e_coeffmat_triangle_q2(points_i)
        call patch_to_k_mat(K,Ke,ids_i)
    end do
end subroutine create_quad_matrix

pure function create_e_coeffmat_triangle_q2(points) result(K_e)
    implicit none
    DOUBLE PRECISION ,INTENT(IN) :: points(3,2)
    DOUBLE PRECISION :: K_e(6,6)
    DOUBLE PRECISION , PARAMETER :: K_before(6,6) = reshape(&
    [1.,1/6.,1/6.,-2/3.,-2/3.,0.,&
    1/6.,0.5,0.,-2/3.,0.,0.,&
    1/6.,0.,0.5,0.,-2/3.,0.,&
    -2/3.,-2/3.,0.,8/3.,0.,-4/3.,&
    -2/3.,0.,-2/3.,0.,8/3.,-4/3.,&
    0.,0.,0.,-4/3.,-4/3.,8/3.]&
    ,[6,6])

    K_e = K_before * jacobian_for_change_variable(points(1:3,:))
end function create_e_coeffmat_triangle_q2

pure function create_fe_vector(points,f) result(f_e)
    implicit none
    DOUBLE PRECISION ,INTENT(IN) :: points(3,2)
    DOUBLE PRECISION , INTENT(IN) :: f(3)
    DOUBLE PRECISION :: f_e(3)
    ! TODO
end function create_fe_vector

pure function jacobian_for_change_variable(points) result(jacobian)
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: points(3,2)
    DOUBLE PRECISION :: jacobian
    jacobian = (points(2,1)-points(1,1))*(points(3,2)-points(1,2)) - (points(3,1)-points(1,1))*(points(2,2)-points(1,2))
end function

subroutine patch_to_k_mat(k,ke,ids)
    ! if ids' size is n, then the shape of Ke should be (n,n)
    implicit none
    DOUBLE PRECISION , INTENT(INOUT) :: K(:,:)
    DOUBLE PRECISION ,INTENT(IN) :: Ke(:,:)
    INTEGER , INTENT(IN) :: ids(:)
    INTEGER i,j
    do j = 1,ubound(ids,1)
        do i = 1,ubound(ids,1)
            K(ids(i),ids(j)) = K(ids(i),ids(j)) + Ke(i,j)
        end do
    end do
end subroutine patch_to_k_mat
end module two_dimentional