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
end module utils
    program main
        ! solve Δu = p 
        !use one_dimensional
        use two_dimentional
        use utils 
        use mkl_service 
        implicit none 
        DOUBLE PRECISION , ALLOCATABLE :: K(:,:),f_vec(:),f_at_points(:)
        type(triangle), ALLOCATABLE :: triangles(:)
        type(node) ,ALLOCATABLE :: nodes(:,:)
        INTEGER n,num_points_in_one_edge, num_points, num_triangles
        Integer , ALLOCATABLE :: dirichlet_ids(:)
        INTEGER , ALLOCATABLE :: ipiv(:)
        INTEGER i,j,count,info
        
        print *,'正方形一辺に何個の三角形の辺を割り当てますか'
        read *,n
        num_points_in_one_edge = 2*n + 1
        num_points = num_points_in_one_edge * num_points_in_one_edge
        num_triangles = 2*n*n
        ALLOCATE(K(num_points,num_points),&
            f_vec(num_points),f_at_points(num_points),&
            triangles(num_triangles),&
            nodes(num_points_in_one_edge,num_points_in_one_edge),&
            dirichlet_ids((num_points_in_one_edge-1)*4),&
            ipiv(num_points))

        ! 正方形の領域にメッシュを作る
        
        count = 1
        do j = 1,num_points_in_one_edge
            do i = 1,num_points_in_one_edge
                nodes(i,j) = node(num_points_in_one_edge*(i-1)+j,[(i-1)*10.d0/2.d0/n,(j-1)*10.d0/2.d0/n])
                if ( i==1 .or. i == num_points_in_one_edge .or. j == 1 .or. j == num_points_in_one_edge ) then
                    dirichlet_ids(count) = num_points_in_one_edge*(i-1)+j
                    count = count + 1
                end if
            end do
        end do
        
        count = 1
        do i = 1,num_points_in_one_edge-2,2
            do j = 1,num_points_in_one_edge-2,2
                triangles(count) = triangle([nodes(i,j),nodes(i,j+2),nodes(i+2,j),nodes(i,j+1),nodes(i+1,j+1),nodes(i+1,j)])
                count = count + 1
            end do
        end do
        do i = 3,num_points_in_one_edge,2
            do j = 3,num_points_in_one_edge,2
                triangles(count) = triangle([nodes(i,j),nodes(i,j-2),nodes(i-2,j),nodes(i,j-1),nodes(i-1,j-1),nodes(i-1,j)])
                count = count + 1
            end do
        end do
        
        f_at_points = 1.d0 ! f=1に固定 ふくざつな関数でもOK

        call create_quad_matrix(K,f_vec,f_at_points,triangles)
        f_vec = - f_vec
        !call print_mat(K)
        do i = 1, ubound(dirichlet_ids,1)
            call dirichlet(K,f_vec,dirichlet_ids(i),0.d0)
        end do
        !call print_mat(K)
        call dgetrf(num_points,num_points,K,num_points,ipiv,info)
        call dgetrs('N',num_points,1, K,num_points,ipiv,f_vec,num_points,info)
        open(18,file = 'ans.csv',status='replace')
        write(18,*) 'X','Y','DATA'
        do i = 1,num_points_in_one_edge
            do j = 1,num_points_in_one_edge
                write(18,*) nodes(i,j)%point(1),',',nodes(i,j)%point(2),',',f_vec(nodes(i,j)%id)
            end do
        end do
        close(18)
        print *,'Hello end'
        ! double precision, allocatable :: a(:,:),b(:,:),f(:)
        ! integer , allocatable :: ipiv(:)
        ! double precision alpha,beta,length
        ! type(line), allocatable :: lines(:)
        ! type(node) :: left,right

        ! alpha = 1.0d0
        ! beta = 1.0d0
        ! read *, n

        ! allocate(a(n+1,n+1),b(n+1,1),f(n+1),lines(n),ipiv(n+1))
        ! f = 1.0d0
        ! length = 1.d0 / n
        ! print *,length
        ! do i = 1,n
        !     left = node(i,length*(i-1),f(i))
        !     right = node(i+1,length * i, f(i+1))
        !     lines(i) = line(left,right)
        ! end do
        ! call make_matrix(A,b(:,1),n,lines)
        ! A = A/dble(n)
        ! b = b/dble(n)
        ! ! ディレクレ
        ! b(1:n+1,1) = b(1:n+1,1) - A(1:n+1,1) * alpha
        ! b(1,1) = alpha
        ! A(1:n+1,1) = 0
        ! A(1,1:n+1) = 0
        ! A(1,1) = 1.0d0

        ! b(1:n+1,1) = b(1:n+1,1) - A(1:n+1,n+1) * 1.d0
        ! b(n+1,1) = 1.0d0
        ! A(1:n+1,n+1) = 0
        ! A(n+1,1:n+1) = 0
        ! A(n+1,n+1) = 1.0d0
        ! !do i = 1, n+1
        ! !    print '(11f7.3)',(A(i,j), j = 1,n+1)
        ! !end do

        ! call dgetrf(n+1,n+1,a,n+1,ipiv,info)
        ! call dgetrs('N',n+1,1, A,n+1,ipiv,b,n+1,info)
        ! open(18,file = 'ans.dat',status='replace')
        ! do i = 1,n+1
        !     write(18,*) b(i,1)
        ! end do
        ! close(18)
    end program
