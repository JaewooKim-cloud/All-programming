subroutine gr(switch,g,nhis,npart,box,ngr,delg,x,rho,ri)

    implicit none
    integer, intent(in) :: nhis,switch,npart
    integer, intent(inout) :: ngr
    real*8, intent(in) :: box, rho
    real*8,dimension(nhis),intent(inout) :: g,ri
    real*8,intent(inout) :: delg
    real*8, dimension(3) :: xr
    real*8, dimension(3,npart), intent(in) :: x
    real*8 :: r,vb,nid
    integer :: i,j,ig
    
    if (switch.eq.0) then
        ngr=0
        delg = box/(2*nhis)
        
        do i=0, nhis
            g(i) = 0
        end do
    else if (switch.eq.1) then
        ngr = ngr + 1
        do i = 1, npart-1
            do j=i+1,npart
                xr(:) = x(:,i)-x(:,j)
                call distpbc3d (xr,box)
                r = sqrt(dot_product(xr,xr))
                if (r.lt.box/2) then
                    ig = int(r/delg)
                    g(ig) = g(ig) + 2
                end if
            end do
        end do
    else if (switch.eq.2) then
        do i=1,nhis
            ri = delg*(i+0.5)
            vb = ((i+1)**3-i**3)*delg**3
            nid = (4.0/3.0)*3.141592*vb*rho
            g(i)=g(i)/(ngr*npart*nid)
        end do
    end if
end subroutine gr
