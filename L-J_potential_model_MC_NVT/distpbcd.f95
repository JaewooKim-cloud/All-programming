subroutine distpbc3d (vec,L)

!Vector should be in the range 0 -> L in all dimensions
!Therefore, we need to apply the following changes if it's not in this range:
    implicit none
    real*8,dimension(3),intent(inout) :: vec
    real*8, intent(in) :: L
    integer :: i
    
    do i=1,3
        if (vec(i).gt.L/2.0) then
            vec(i) = vec(i)-L
        end if
        
        if (vec(i).lt.-L/2.0) then
            vec(i) = vec(i)+L
        end if
    
    end do
    





end subroutine distpbc3d
