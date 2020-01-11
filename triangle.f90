module two_dimentional
implicit none
type node
    INTEGER id
    DOUBLE PRECISION point(2)
end type
type triangle
    TYPE(node) :: nodes(6)
end type

contains
subroutine create_quad_matrix(K_mat,f_vec,f_at_points,triangles)
    implicit none
    !INTEGER , INTENT(IN) :: n ! size of matrix A(n*n) and vector b
    DOUBLE PRECISION , INTENT(OUT) :: K_mat(:,:),f_vec(:)
    DOUBLE PRECISION , INTENT(IN) :: f_at_points(:)
    TYPE(triangle) :: triangles(:),triangle_i
    DOUBLE PRECISION :: points_i(6,2),Ke(6,6),fe(3)
    INTEGER :: ids_i(6)
    INTEGER i,l

    K_mat(:,:) = 0
    f_vec(:) = 0

    do i = 1,ubound(triangles,1)
        triangle_i = triangles(i)
        ids_i = triangle_i%nodes%id
        do l = 1,6
            points_i(l,:) = triangle_i%nodes(l)%point
        end do
        Ke = create_e_coeffmat_triangle_q2(points_i)
        fe = create_fe_vector(points_i(1:3,:),f_at_points(ids_i(1:3)))
        call patch_to_k_mat(K_mat,Ke,ids_i)
        call patch_to_f_vec(f_vec,fe,ids_i(1:3))
    end do
end subroutine create_quad_matrix

pure function create_e_coeffmat_triangle_q2(points) result(K_e)
    implicit none
    DOUBLE PRECISION ,INTENT(IN) :: points(3,2)
    DOUBLE PRECISION :: K_e(6,6)
    DOUBLE PRECISION , PARAMETER :: phi_xi_phi_xi(6,6) = reshape(&
    [0.5, 0.166666666666667, 0., -0.666666666666667, 0., 0.,&
       0.166666666666667, 0.5, 0., -0.666666666666667, 0., 0.,&
       0., 0., 0., 0., 0., 0., -0.666666666666667, -0.666666666666667, 0.,&
       1.33333333333333, 0., 0., 0., 0., 0., 0., 1.33333333333333,&
       -1.33333333333333, 0., 0., 0., 0., -1.33333333333333, 1.33333333333333],[6,6])
    DOUBLE PRECISION , PARAMETER :: phi_xi_phi_eta(6,6) = reshape(&
    [1.00000000000000, 0.166666666666667, 0.166666666666667,&
       -0.666666666666667, 0., -0.666666666666667, 0.166666666666667, 0.,&
       -0.166666666666667, -0.666666666666667, 0.666666666666667, 0.,&
       0.166666666666667, -0.166666666666667, 0., 0., 0.666666666666667,&
       -0.666666666666667, -0.666666666666667, -0.666666666666667, 0.,&
       1.33333333333333, -1.33333333333333, 1.33333333333333, 0.,&
       0.666666666666667, 0.666666666666667, -1.33333333333333,&
       1.33333333333333, -1.33333333333333, -0.666666666666667, 0.,&
       -0.666666666666667, 1.33333333333333, -1.33333333333333,&
       1.33333333333333],[6,6])
    DOUBLE PRECISION, PARAMETER :: phi_eta_phi_eta(6,6) = reshape(&
        [0.5, 0., 0.166666666666667, 0., 0., -0.666666666666667,&
       0., 0., 0., 0., 0., 0., 0.166666666666667, 0., 0.500000000000000, 0., 0.,&
       -0.666666666666667, 0., 0., 0., 1.33333333333333, -1.33333333333333,&
       0., 0., 0., 0., -1.33333333333333, 1.33333333333333, 0.,&
       -0.666666666666667, 0., -0.666666666666667, 0., 0., 1.33333333333333],&
    [6,6])
    DOUBLE PRECISION x2(2),x3(2)
    x2 = points(2,:)-points(1,:)
    x3 = points(3,:)-points(1,:)
    K_e = (phi_xi_phi_xi*DOT_PRODUCT(x2,x2)-phi_xi_phi_eta*DOT_PRODUCT(x2,x3)+phi_eta_phi_eta*DOT_PRODUCT(x3,x3)) &
    /jacobian_for_change_variable(points(1:3,:))
end function create_e_coeffmat_triangle_q2

pure function create_fe_vector(points,f) result(f_e)
    implicit none
    DOUBLE PRECISION ,INTENT(IN) :: points(3,2)
    DOUBLE PRECISION , INTENT(IN) :: f(3)
    DOUBLE PRECISION :: f_e(3)
    DOUBLE PRECISION ,PARAMETER :: mat(3,3) = reshape(&
    ![11/12.,-0.25,-0.875,-0.25,0.125,5./24,-0.875,5./24,11./12]&
    [1/12.,1/24.,1/24.,1/24.,1/12.,1/24.,1/24.,1/24.,1/12.]&
    ,[3,3])
    f_e = matmul(mat,f)*jacobian_for_change_variable(points)!*jacobian_for_change_variable(points)

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
subroutine patch_to_f_vec(f_vec,fe,ids)
    ! if ids' size is n, then the shape of Ke should be (n,n)
    implicit none
    DOUBLE PRECISION , INTENT(INOUT) :: f_vec(:)
    DOUBLE PRECISION ,INTENT(IN) :: fe(:)
    INTEGER , INTENT(IN) :: ids(:)
    INTEGER i
    do i = 1,3
        f_vec(ids(i)) = f_vec(ids(i)) + fe(i)
    end do
end subroutine patch_to_f_vec
end module two_dimentional