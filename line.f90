! ∇^2u = f
module one_dimensional
    implicit none
    type node
        integer id
        double precision point
        double precision f
    end type

    type line
        type(node) left
        type(node) right
    end type

    contains
        subroutine make_matrix(A,b,n,lines)
            implicit none
            integer,intent(in) :: n
            type(line),intent(in) :: lines(n)
            double precision,intent(out) :: A(n+1,n+1),b(n+1)
            integer i
            integer left_id,right_id
            double precision distance,reciprocal_distance, left_f, right_f
            type(line) :: line_i

            A = 0.d0
            b = 0.d0
            do i = 1,n
                line_i = lines(i)
                left_id = line_i%left%id
                right_id = line_i%right%id
                left_f = line_i%left%f
                right_f = line_i%right%f
                distance = line_i%right%point - line_i%left%point
                reciprocal_distance = 1.d0/ distance
                ! 合成行列に代入
                A(left_id,left_id) = A(left_id,left_id) + reciprocal_distance
                A(left_id,right_id) = A(left_id,right_id) - reciprocal_distance
                A(right_id,left_id) = A(right_id,left_id) - reciprocal_distance
                A(right_id,right_id) = A(right_id,right_id) + reciprocal_distance
                b(left_id) = b(left_id) - (2.d0 * left_f + right_f)*distance/6.0d0
                b(right_id) = b(right_id) - (left_f + 2.d0 * right_f)*distance / 6.0d0
            end do
        end subroutine make_matrix

end module one_dimensional