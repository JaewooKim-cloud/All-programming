program md_basic

          
        
    implicit none
    integer,parameter :: seed = 864412
    integer, parameter :: npart = 125, nsteps = 10000, nhis = 1000
    real*8 :: density , temp, dt, mass,tmax,t
    real*8, dimension (3,npart) :: x, xm,f,v
    real*8, dimension (nhis) :: ri,g
    real*8 :: L, etot,en,rc2,ecut,rc,delg
    integer :: coun,ngr,i
    
    coun = 0
    density = 0.85
    temp = 2.0
    dt = 0.001
    tmax=10
    mass = 1.0
    
    call srand(seed)
    call initcubicgrid(npart,density,L,x) ! place the particle on a lattice
    call init (npart,density,temp,L,x,xm,dt,v) ! set initial configuraion
    rc = 0.5*L
    rc2 = (rc)**2
    ecut = 4*(1.0/(rc**12)-1.0/(rc**6))
    
    call gr(0,g,nhis,npart,L,ngr,delg,x,density,ri)
   t = 0
   
   !start MD loop
    do while (t.lt.tmax)
        
        call force(x,f,en,npart,L,rc2,ecut)
        call integrate(f,en,x,xm,etot,dt,temp,npart,v)
        t = t + dt
        
       
      if (coun.gt.1000) then
           call gr(1,g,nhis,npart,L,ngr,delg,x,density,ri)
    end if
       
    
       coun = coun+1
    end do
    
    call gr(2,g,nhis,npart,L,ngr,delg,x,density,ri)
     open(unit=1,file="rad.txt")
     10 format(f10.4, 5x , f10.4)
    do i=1,nhis
        write(unit=1,fmt=10),ri(i), g(i)
    end do
    close(unit=1)
    
end program md_basic


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

subroutine integrate(f,en,x,xm,etot,delt,temp,npart,v)

    implicit none
    real*8, dimension(3,npart),intent(inout) :: x,f,xm,v
    real*8, intent(in) :: en,delt
    real*8, dimension(3) :: xx,sumv
    integer, intent(in) :: npart
    real*8, intent(inout) :: etot,temp
    real*8 :: sumv2
    integer :: i
    
    
    sumv = 0
    sumv2 = 0
    
    do i=1,npart
        xx(:) = 2*x(:,i) - xm(:,i) + (delt**2)*f(:,i) !verlet algorithm
        v(:,i) = (xx(:)-xm(:,i))/(2*delt)
        sumv(:) = sumv(:) + v(:,i)
        sumv2 = sumv2 + dot_product(v(:,i),v(:,i)) ! total kinetic energy
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
            if (r2.lt.rc2) then !test cut-of
                r2i = 1.0/r2
                r6i = r2i**3
                ff = 48*r2i*r6i*(r6i-0.5) !LJ force
                f(:,i) = f(:,i) + ff*xr(:)
                f(:,j) = f(:,j) - ff*xr(:) 
                en = en + 4*r6i*(r6i-1) - ecut
            end if           
        end do
    end do 
end subroutine force

subroutine init (npart,density,temp,L,x,xm,dt,v)
    
    implicit none
    integer, intent(in) :: npart
    real*8, intent(in) :: density, temp,dt
    real*8, intent(inout) :: L
    real*8, dimension(3,npart), intent(inout) :: x,xm,v
    real*8,dimension(3) :: sumv
    real*8 :: sumv2,fs
    integer :: i
    
    sumv=0
    sumv2=0
    
    
    
    do i=1, npart
        v(1,i) = (rand()-0.5) ! give random velocities
        v(2,i) = (rand()-0.5) ! give random velocities
        v(3,i) = (rand()-0.5) ! give random velocities
         
        sumv(:) = sumv(:) + v(:,i) ! velocity center of mass
        sumv2 = sumv2 + dot_product(v(:,i),v(:,i))    !kinetic energy
        
    end do
    
    sumv(:) = sumv(:)/npart
    
    sumv2 = sumv2/npart !mean squared velovity
    fs = sqrt(3*temp/sumv2)!scale factor of velocity
     
    do i=1, npart
        v(:,i)=(v(:,i)-sumv(:))*fs ! velocity center of mass to zero
        xm(:,i) = x(:,i)-v(:,i)*dt
    end do
  
    call pbc3d (xm,L,npart) !periodic bdry condition
    
end subroutine init

subroutine pbc3d (vec,L,npart)

!Vector should be in the range 0 -> L in all dimensions
!Therefore, we need to apply the following changes if it's not in this range:
    implicit none
    integer, intent(in) :: npart
    real*8,dimension(3,npart),intent(inout) :: vec
    real*8, intent(in) :: L
    integer :: i,j
    
    
    do j=1,npart
        do i=1,3
            if (vec(i,j).gt.L) then 
              
                vec(i,j) = vec(i,j)-L
            end if
       
            if (vec(i,j).lt.0) then 
                vec(i,j) = vec(i,j)+L
            end if
             
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
   
        if (vec(i).gt.L/2.0) then
            vec(i) = vec(i)-L
            
        end if
         
        if (vec(i).lt.-L/2.0) then
       
            vec(i) = vec(i)+L
        end if
         
    end do
    
end subroutine distpbc3d

subroutine gr(switch,g,nhis,npart,L,ngr,delg,x,density,ri)

    implicit none
    integer, intent(in) :: nhis,switch,npart
    integer, intent(inout) :: ngr
    real*8, intent(in) :: L, density
    real*8,dimension(nhis),intent(inout) :: g,ri
    real*8,intent(inout) :: delg
    real*8, dimension(3) :: xr
    real*8, dimension(3,npart), intent(in) :: x
    real*8 :: r,vb,nid
    integer :: i,j,ig
    
    if (switch.eq.0) then
        ngr=0
        delg = L/(2*nhis)
        
        do i=1, nhis
            g(i) = 0
        end do
    else if (switch.eq.1) then
        ngr = ngr + 1
        do i = 1, npart-1
            do j=i+1,npart
                xr(:) = x(:,i)-x(:,j)
                call distpbc3d (xr,L)
                r = (dot_product(xr,xr))**0.5
          
                if (r.lt.L/2.0) then
                    ig = int(r/delg)
                    g(ig) = g(ig) + 2
                  
                end if
            end do
        end do
    else if (switch.eq.2) then
        do i=1,nhis
            ri(i) = delg*(i+0.5)
            vb = ((i+1)**3-i**3)*delg**3
            nid = (4.0/3.0)*3.141592*vb*density
            g(i)=g(i)/(ngr*npart*nid)
        end do
    end if
end subroutine gr

