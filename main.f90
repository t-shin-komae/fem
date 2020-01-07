    program main
        ! solve Δu = p where (0,1) u(0)= alpha, du/dx |x=1 = beta
        use one_dimensional
        use mkl_service
        implicit none
        integer i,j,n,info
        double precision, allocatable :: a(:,:),b(:,:),f(:)
        integer , allocatable :: ipiv(:)
        double precision alpha,beta,length
        type(line), allocatable :: lines(:)
        type(node) :: left,right

        alpha = 1.0d0
        beta = 1.0d0
        read *, n

        allocate(a(n+1,n+1),b(n+1,1),f(n+1),lines(n),ipiv(n+1))
        f = 1.0d0
        length = 1.d0 / n
        print *,length
        do i = 1,n
            left = node(i,length*(i-1),f(i))
            right = node(i+1,length * i, f(i+1))
            lines(i) = line(left,right)
        end do
        call make_matrix(A,b(:,1),n,lines)
        A = A/dble(n)
        b = b/dble(n)
        ! ディレクレ
        b(1:n+1,1) = b(1:n+1,1) - A(1:n+1,1) * alpha
        b(1,1) = alpha
        A(1:n+1,1) = 0
        A(1,1:n+1) = 0
        A(1,1) = 1.0d0

        b(1:n+1,1) = b(1:n+1,1) - A(1:n+1,n+1) * 1.d0
        b(n+1,1) = 1.0d0
        A(1:n+1,n+1) = 0
        A(n+1,1:n+1) = 0
        A(n+1,n+1) = 1.0d0
        !do i = 1, n+1
        !    print '(11f7.3)',(A(i,j), j = 1,n+1)
        !end do

        call dgetrf(n+1,n+1,a,n+1,ipiv,info)
        call dgetrs('N',n+1,1, A,n+1,ipiv,b,n+1,info)
        open(18,file = 'ans.dat',status='replace')
        do i = 1,n+1
            write(18,*) b(i,1)
        end do
        close(18)
    end program