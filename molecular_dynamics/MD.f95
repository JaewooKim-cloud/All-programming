module MD

contains

subroutine initcubicgrid(npart,density,L,coords)

    implicit none
    real*8, dimension (3,npart), intent(inout) :: coords
    real*8, intent(in) :: density
    real*8, intent(inout) :: L
    integer :: i
    integer, intent(in) :: npart
    real, dimension(3) :: in_dex
    real :: ncube
    
    ncube = 0
    
    do while (ncube**(3).lt.npart)
        ncube = ncube+1
    end do
    
    L = (npart/density)**(1.0/3.0)
    
    coords(:,:) = 0
    
    in_dex(:)=0
    
    do i=1,npart
        coords(1,i) = (in_dex(1)+0.5)*(L/ncube)
        coords(2,i) = (in_dex(2)+0.5)*(L/ncube)
        coords(3,i) = (in_dex(3)+0.5)*(L/ncube)
        
        in_dex(1)=in_dex(1)+1
        if(in_dex(1).eq.ncube) then
            in_dex(1)=0
            in_dex(2)=in_dex(2)+1
            if(in_dex(2).eq.ncube) then
                in_dex(2)=0
                in_dex(3)=in_dex(3)+1
            end if
        end if
    end do
    
end subroutine initcubicgrid

subroutine integrate(f,en,x,xm,etot,delt,temp,npart)

    implicit none
    real*8, dimension(3,npart),intent(inout) :: x,f,xm
    real*8, intent(in) :: en,delt
    real*8, dimension(3) :: xx,vi,sumv
    integer, intent(in) :: npart
    real*8, intent(inout) :: etot,temp
    real*8 :: sumv2
    integer :: i
    
    
    sumv = 0
    sumv2 = 0
    
    do i=1,npart
        xx(:) = 2*x(:,i) - xm(:,i) + (delt**2)*f(:,i) !verlet algorithm
        vi(:) = (xx(:)-xm(:,i))/(2*delt)
        sumv(:) = sumv(:) + vi(:)
        sumv2 = sumv2 + dot_product(vi,vi) ! total kinetic energy
        xm(:,i) = x(:,i) !update positions previous time
        x(:,i) = xx(:) !update positions current time 
    end do
    
    temp = sumv2/(3.0*npart) !instantaneous temperature
    
    etot = (en + 0.5*sumv2)/npart

end subroutine integrate

subroutine force(x,f,en,npart,box,rc2,ecut)

    implicit none
    integer :: i,j
    integer, intent(in) :: npart
    real*8, dimension(3) :: xr
    real*8, dimension(3,npart),intent(inout) :: x,f
    real*8, intent(in) :: box
    real*8 :: r2,ff,r2i,r6i
    real*8, intent(in) :: rc2,ecut
    real*8, intent(inout) :: en
    
    en = 0
    f(:,:) = 0 !set forces to zero
    
    do i=1, npart-1
        do j=i+1, npart
            xr(:) = x(:,i)-x(:,j) !xr is diffrence of positions of i and j particle
            call distpbc3d (xr,box) ! periodic bdry conditions
            r2 = dot_product(xr,xr)
            if (r2.lt.rc2) then !test cut-off
                r2i = 1.0/r2
                r6i = r2i**3
                ff = 48*r2i*r6i*(r6i-0.5) !LJ force
                f(:,i) = f(:,i) + ff*xr(:)
                f(:,j) = f(:,j) + ff*xr(:) 
                en = en + 4*r6i*(r6i-1) - ecut
            end if            
        end do
    end do 
end subroutine force

subroutine init (npart,density,temp,L,x,xm,dt)
    
    implicit none
    integer, intent(in) :: npart
    real*8, intent(in) :: density, temp,dt
    real*8, intent(inout) :: L
    real*8, dimension(3,npart), intent(inout) :: x,xm
    real*8,dimension(3) :: sumv,v
    real*8 :: sumv2,fs
    integer :: i
    
    sumv=0
    sumv2=0
    
    call initcubicgrid(npart,density,L,x) ! place the particle on a lattice
    
    do i=1, npart
        v(1) = (rand()-0.5) ! give random velocities
        v(2) = (rand()-0.5) ! give random velocities
        v(3) = (rand()-0.5) ! give random velocities
        
        sumv(:) = sumv(:) + v(:) ! velocity center of mass
        sumv2 = sumv2 + v(1)**2    !kinetic energy
        sumv2 = sumv2 + v(2)**2    !kinetic energy
        sumv2 = sumv2 + v(3)**2    !kinetic energy
    end do
    sumv(:) = sumv(:)/npart
    sumv2 = sumv2/npart !mean squared velovity
    fs = sqrt(3*temp/sumv2)!scale factor of velocity
    do i=1, npart
        v(:)=(v(:)-sumv(:))*fs ! velocity center of mass to zero
        xm(:,i) = x(:,i)-v(:)*dt
    end do
    call pbc3d (xm,L,npart) !periodic bdry condition
    
end subroutine init

subroutine pbc3d (vec,L,npart)

!Vector should be in the range 0 -> L in all dimensions
!Therefore, we need to apply the following changes if it's not in this range:
    implicit none
    real*8,dimension(3,npart),intent(inout) :: vec
    real*8, intent(in) :: L
    integer :: i,j
    integer, intent(in) :: npart
    
    do j=1,npart
        do i=1,3
            do while (vec(i,j).gt.L) 
                vec(i,j) = vec(i,j)-L
            end do
        
            do while (vec(i,j).lt.0) 
                vec(i,j) = vec(i,j)+L
            end do
        end do
    end do
    
end subroutine pbc3d



subroutine distpbc3d (vec,L)

!Vector should be in the range -L/2-> L/2 in all dimensions
!Therefore, we need to apply the following changes if it's not in this range:
    implicit none
    real*8,dimension(3),intent(inout) :: vec
    real*8, intent(in) :: L
    integer :: i
    
    do i=1,3
        do while (vec(i).gt.L/2.0) 
            vec(i) = vec(i)-L
        end do
        
        do while (vec(i).lt.-L/2.0)
            vec(i) = vec(i)+L
        end do
    end do
    
end subroutine distpbc3d

end module MD
