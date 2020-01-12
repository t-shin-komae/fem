module two_dimentional
implicit none
type node
    INTEGER id
    DOUBLE PRECISION point(2)
end type
type triangle
    TYPE(node) :: nodes(6)
end type
type node_with_val
    INTEGER id
    DOUBLE PRECISION point(2)
    DOUBLE PRECISION u
end type
type triangle_with_val
    TYPE(node_with_val) :: nodes(6)
end type

contains
subroutine create_quad_matrix(K_mat,f_vec,f_at_points,triangles)
    implicit none
    !INTEGER , INTENT(IN) :: n ! size of matrix A(n*n) and vector b
    DOUBLE PRECISION , INTENT(OUT) :: K_mat(:,:),f_vec(:)
    DOUBLE PRECISION , INTENT(IN) :: f_at_points(:)
    TYPE(triangle) :: triangles(:),triangle_i
    DOUBLE PRECISION :: points_i(6,2),Ke(6,6),fe(6)
    INTEGER :: ids_i(6)
    INTEGER i,l

    K_mat(:,:) = 0
    f_vec(:) = 0

    do i = 1,ubound(triangles,1)
        triangle_i = triangles(i)
        ids_i = triangle_i%nodes%id
        do l = 1,6
            points_i(l,:) = triangle_i%nodes(l)%point
            ids_i(l) = triangle_i%nodes(l)%id
        end do
        Ke = create_e_coeffmat_triangle_q2(points_i(1:3,:))
        ! fe = create_fe_vector(points_i(1:3,:),[f_at_points(ids_i(1)),f_at_points(ids_i(2)),f_at_points(ids_i(3))])
        fe = create_fe_vector_q2(points_i,f_at_points(ids_i))
        call patch_to_k_mat(K_mat,Ke,ids_i)
        ! call patch_to_f_vec(f_vec,fe,ids_i(1:3))
        call patch_to_f_vec_q2(f_vec,fe,ids_i)
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
    K_e = (phi_xi_phi_xi*DOT_PRODUCT(x3,x3)-phi_xi_phi_eta*DOT_PRODUCT(x2,x3)+phi_eta_phi_eta*DOT_PRODUCT(x2,x2)) &
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
    f_e = matmul(mat,f)*jacobian_for_change_variable(points)

end function create_fe_vector

pure function create_fe_vector_q2(points,f) result(f_e)
    implicit none
    DOUBLE PRECISION ,INTENT(IN) :: points(6,2)
    DOUBLE PRECISION , INTENT(IN) :: f(6)
    DOUBLE PRECISION :: f_e(6)
    DOUBLE PRECISION ,PARAMETER :: mat(6,6) = reshape(&
    [2.41666666666667, 0.597222222222222, 0.597222222222222,&
       -2.00000000000000, 0.388888888888889, -2.00000000000000,&
       0.597222222222222, 0.216666666666667, 0.0638888888888889,&
       -0.666666666666667, 0.0666666666666667, -0.277777777777778,&
       0.597222222222222, 0.0638888888888889, 0.216666666666667,&
       -0.277777777777778, 0.0666666666666667, -0.666666666666667,&
       -2.00000000000000, -0.666666666666667, -0.277777777777778,&
       2.22222222222222, -0.222222222222222, 1.11111111111111,&
       0.388888888888889, 0.0666666666666667, 0.0666666666666667,&
       -0.222222222222222, 0.0888888888888889, -0.222222222222222,&
       -2.00000000000000, -0.277777777777778, -0.666666666666667,&
       1.11111111111111, -0.222222222222222, 2.22222222222222],[6,6])
    f_e = matmul(mat,f)*jacobian_for_change_variable(points(1:3,:))

end function create_fe_vector_q2

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

subroutine patch_to_f_vec_q2(f_vec,fe,ids)
    ! if ids' size is n, then the shape of Ke should be (n,n)
    implicit none
    DOUBLE PRECISION , INTENT(INOUT) :: f_vec(:)
    DOUBLE PRECISION ,INTENT(IN) :: fe(:)
    INTEGER , INTENT(IN) :: ids(:)
    INTEGER i
    do i = 1,6
        f_vec(ids(i)) = f_vec(ids(i)) + fe(i)
    end do
end subroutine patch_to_f_vec_q2

pure elemental function is_include(tri,x,y)
    implicit none
    TYPE(triangle_with_val) , INTENT(IN) :: tri
    DOUBLE PRECISION , INTENT(IN) :: x,y
    LOGICAL is_include
    DOUBLE PRECISION p1(2),p2(2),p3(2),point(2),coeff(2),s,t
    p1 = tri%nodes(1)%point
    p2 = tri%nodes(2)%point
    p3 = tri%nodes(3)%point
    point = [x,y] - p1
    p2 = p2-p1
    p3 = p3-p1
    coeff = MATMUL(inverse_2([p2,p3]),point)
    s = coeff(1)
    t = coeff(2)
    if ( 0.<= s .and. s <= 1. .and. 0. <= t .and. t <= 1. .and. 0.<=s+t .and. s+t<1.) then
        is_include = .true.
    else
        is_include = .false.
    end if
end function is_include
pure function inverse_2(mat)
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: mat(2,2)
    DOUBLE PRECISION inverse_2(2,2)
    inverse_2(1,1) = mat(2,2)
    inverse_2(1,2) = -mat(1,2)
    inverse_2(2,1) = -mat(2,1)
    inverse_2(2,2) = mat(1,1)
    inverse_2 = inverse_2/(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))
end function inverse_2

pure function calc_val_in_triangle(tri,x,y)
    implicit none
    TYPE(triangle_with_val) ,INTENT(IN) :: tri
    DOUBLE PRECISION ,INTENT(IN) :: x,y
    DOUBLE PRECISION calc_val_in_triangle
    DOUBLE PRECISION xi(2)
    DOUBLE PRECISION p1(2),p2(2),p3(2),point(2)
    p1 = tri%nodes(1)%point
    p2 = tri%nodes(2)%point
    p3 = tri%nodes(3)%point
    point = [x,y] - p1
    p2 = p2-p1
    p3 = p3-p1
    xi= MATMUL(inverse_2([p2,p3]),point)
    
    calc_val_in_triangle =&
        tri%nodes(1)%u*phi1(xi(1),xi(2))+&
        tri%nodes(2)%u*phi2(xi(1),xi(2))+&
        tri%nodes(3)%u*phi3(xi(1),xi(2))+&
        tri%nodes(4)%u*phi4(xi(1),xi(2))+&
        tri%nodes(5)%u*phi5(xi(1),xi(2))+&
        tri%nodes(6)%u*phi6(xi(1),xi(2))
    
end function

pure function phi1(xi,eta) 
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: xi,eta
    DOUBLE PRECISION phi1
    phi1 = xi*(2.*xi+4.*eta-3.) + eta*(2.*eta-3.) +1.
end function
pure function phi2(xi,eta) 
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: xi,eta
    DOUBLE PRECISION phi2 
    phi2 = xi*(2*xi-1)
end function
pure function phi3(xi,eta) 
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: xi,eta
    DOUBLE PRECISION phi3
    phi3 = eta*(2*eta-1)
end function
pure function phi4(xi,eta) 
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: xi,eta
    DOUBLE PRECISION phi4
    phi4 = 4*(1-xi-eta)*xi
end function
pure function phi5(xi,eta) 
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: xi,eta
    DOUBLE PRECISION phi5
    phi5 = 4*xi*eta
end function
pure function phi6(xi,eta) 
    implicit none
    DOUBLE PRECISION , INTENT(IN) :: xi,eta
    DOUBLE PRECISION phi6
    phi6 = 4*(1-xi-eta)*eta
end function
pure function create_triangle_with_val(tri,u) result(tri_with_val)
    implicit none
    TYPE(triangle),INTENT(IN) :: tri
    DOUBLE PRECISION,INTENT(IN) :: u(:)
    TYPE(triangle_with_val) tri_with_val
    TYPE(node_with_val) nodes(6)
    INTEGER i
    do i = 1,6
        nodes(i) = create_node_with_val(tri%nodes(i),u)
    end do
    tri_with_val%nodes = nodes
end function

pure function create_node_with_val(node_i,u) result(node_with_val_i)
    implicit none
    TYPE(node) , INTENT(IN) :: node_i
    DOUBLE PRECISION , INTENT(IN) :: u(:)
    TYPE(node_with_val) node_with_val_i
    node_with_val_i%point = node_i%point
    node_with_val_i%id = node_i%id
    node_with_val_i%u = u(node_i%id)
end function

end module two_dimentional