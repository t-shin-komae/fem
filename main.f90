    program main
        ! solve Δu = p 
        !use one_dimensional
        use two_dimentional
        use utils 
        use mkl_service 
        use mesh
        implicit none 
        integer , parameter :: resolution = 64
        DOUBLE PRECISION , ALLOCATABLE :: K(:,:),f_vec(:),f_at_points(:),dirichlet_val(:)
        type(triangle), ALLOCATABLE :: triangles(:)
        type(node) ,ALLOCATABLE :: nodes(:,:)
        INTEGER n,num_points_in_one_edge, num_points, num_triangles
        Integer , ALLOCATABLE :: dirichlet_ids(:)
        INTEGER , ALLOCATABLE :: ipiv(:)
        INTEGER i,j,l,count,info
        DOUBLE PRECISION u(resolution+1,resolution+1)
        type(triangle_with_val),allocatable :: tri_with_val(:)
        
        print *,'正方形一辺に何個の三角形の辺を割り当てますか'
        read *,n
        num_points_in_one_edge = 2*n + 1
        num_points = num_points_in_one_edge * num_points_in_one_edge
        num_triangles = 2*n*n
        ALLOCATE(K(num_points,num_points),&
            f_vec(num_points),f_at_points(num_points),&
            triangles(num_triangles),&
            tri_with_val(num_triangles),&
            nodes(num_points_in_one_edge,num_points_in_one_edge),&
            dirichlet_ids((num_points_in_one_edge-1)*4),&
            dirichlet_val((num_points_in_one_edge-1)*4),&
            ipiv(num_points))

        ! 正方形の領域にメッシュを作る
        
        count = 1
        do j = 1,num_points_in_one_edge
            do i = 1,num_points_in_one_edge
                nodes(i,j) = node(id = num_points_in_one_edge*(i-1)+j,point = [(i-1)*1.d0/2.d0/n,(j-1)*1.d0/2.d0/n])
                if (j == 1 ) then
                    dirichlet_ids(count) = num_points_in_one_edge*(i-1)+j
                    dirichlet_val(count) = 0.d0
                    count = count + 1
                else if (j == num_points_in_one_edge) then
                    dirichlet_ids(count) = num_points_in_one_edge*(i-1)+j
                    dirichlet_val(count) = 0.d0
                    count = count + 1
                else if ( i==1 .or. i == num_points_in_one_edge  ) then
                    dirichlet_ids(count) = num_points_in_one_edge*(i-1)+j
                    dirichlet_val(count) = 0.d0
                    count = count + 1
                end if
            end do
        end do
        
        count = 1
        do i = 1,num_points_in_one_edge-2,2
            do j = 1,num_points_in_one_edge-2,2
                triangles(count) = triangle([nodes(i,j),nodes(i+2,j),nodes(i,j+2),nodes(i+1,j),nodes(i+1,j+1),nodes(i,j+1)])
                count = count + 1
            end do
        end do
        do i = 3,num_points_in_one_edge,2
            do j = 3,num_points_in_one_edge,2
                triangles(count) = triangle([nodes(i,j),nodes(i-2,j),nodes(i,j-2),nodes(i-1,j),nodes(i-1,j-1),nodes(i,j-1)])
                count = count + 1
            end do
        end do
        
        f_at_points = -1.d0 ! f=1に固定 ふくざつな関数でもOK

        call create_quad_matrix(K,f_vec,f_at_points,triangles)
        f_vec = - f_vec
        !call print_mat(K)
        do i = 1, ubound(dirichlet_ids,1)
            call dirichlet(K,f_vec,dirichlet_ids(i),dirichlet_val(i))
        end do
        call dgetrf(num_points,num_points,K,num_points,ipiv,info)
        call dgetrs('N',num_points,1, K,num_points,ipiv,f_vec,num_points,info)

        do i = 1, ubound(triangles,1)
            tri_with_val(i) = create_triangle_with_val(triangles(i),f_vec)
        end do
        open(18,file = 'ans.csv',status='replace')
        write(18,*) 'X',',','Y',',','DATA'
        do j = 0,resolution
            do i = 0,resolution
                do l = 1,ubound(triangles,1)
                    if (is_include(tri_with_val(l),i/dble(resolution),j/dble(resolution))) then
                        write(18,*) i/dble(resolution),',',j/dble(resolution),',',calc_val_in_triangle(tri_with_val(l),i/dble(resolution),j/dble(resolution))!,',',triangles(l)%nodes(1:3)
                        exit
                    end if
                end do
            end do
            write(18,*) ' '
        end do
        close(18)
        print *,'Hello end'
    end program
