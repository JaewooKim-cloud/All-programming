program mc_micro

          
            
    implicit none
    integer,parameter :: seed = 12451
    integer, parameter ::  nsteps = 30000, printFreq = 100
    real*8 :: density ,temp, beta,maxDr,zz
    integer :: i,j
    integer :: npart = 108
    real*8, dimension (3,300) :: coords
 
    real*8 :: L,deltaE, pid, energy,V
    real*8,dimension(3) :: rtrial
    real*8,external :: LJ_Energy, vir
    real*8,external :: LJ_EnergyChange
    real*8, dimension(20000) :: density1
   
 
    V = 180.0
    maxDr = 0.5
    pid = 3.0
    temp = 2.0
    beta = 1.0/temp
   
    zz = beta * pid
    
    call srand(seed)
    
   call initcubicgrid(npart,density,L,coords,V) ! set initial configuraion
   energy = LJ_Energy(coords,L,npart,density)
    
  ! zz is defined as         zz   = exp(beta*mu^b)/Lamda^3
      !                                = beta*pid
!        excess chemical pot.     muex = mu^b -mu^id
!                                     = mu^0 + ln(beta*pid)/beta - mu^0 - ln(rho)
!                                      = ln(zz)/beta - ln <rho>

!pid     = ideal gas pressure reservoir (transfered to chemical potential)
! in this simulation, we set T=2.0, N=108 for initial N and V = 180.0
   ! mc cycle
    do i=1, nsteps
        do j=1, npart
            
                rtrial(1)=coords(1,j)+maxDr*(rand()-0.5)
                rtrial(2)=coords(2,j)+maxDr*(rand()-0.5)
                rtrial(3)=coords(3,j)+maxDr*(rand()-0.5)
               call pbc3d (rtrial,L) ! apply periodic boundary condition
                deltaE = LJ_EnergyChange(coords,rtrial,j,L,npart)
                if(rand().lt.exp(-beta*deltaE)) then
                    coords(1,j) = rtrial(1)
                    coords(2,j) = rtrial(2)
                    coords(3,j) = rtrial(3)                    
                end if
        end do
        call mcexc (L,V,npart,coords,zz,beta) ! for every one step, we perform mcexc
        if (i.gt.10000) then
            density1(i-10000) = npart/(L**3)
        end if
    end do
    density = 0
    
   
    open(unit=1,file="density1.txt")
     10 format(i5, 5x , f10.4,5x)
    do i=1,20000
        write(unit=1,fmt=10),i, density1(i)
    end do
    close(unit=1)
    
    do i=1,20000
        density = density + density1(i)
    end do
    
    density = density/20000.0
    
    print *, density
   
  
 
    
  
  
 
end program mc_micro

real*8 function vir(npart,coords,L)

    implicit none
    integer, intent(in) :: npart
    real*8, dimension(3,300) ,intent(in) :: coords
    integer i,j
    real*8,dimension(3) ::r_ij
    real*8, intent(in) :: L
    
    vir = 0
    
    do i=1,npart-1
        do j=i+1,npart
            
            r_ij(1)=coords(1,i)-coords(1,j)
            r_ij(2)=coords(2,i)-coords(2,j)
            r_ij(3)=coords(3,i)-coords(3,j)
            
            call distpbc3d (r_ij,L)
           
           if((dot_product(r_ij,r_ij))**(0.5).lt.0.5*L) then
           vir = vir + 48*(dot_product(r_ij,r_ij)**(-6.0)-0.5*(dot_product(r_ij,r_ij)**(-3.0))) 
            end if
        end do
    
    end do
    vir = vir*(1.0/3.0)
end function vir

subroutine initcubicgrid(npart,density,L,coords,V)

    implicit none
    real*8, dimension (3,300), intent(inout) :: coords
    real*8, intent(in) :: V
    real*8, intent(inout) :: L,density
    integer :: i
    integer, intent(in) :: npart
    real, dimension(3) :: in_dex
    real*8 :: ncube
    
    ncube = 0
    
    do while (ncube**(3).lt.npart)
        ncube = ncube+1
    end do
    
    density = npart/V
   
    L = (V)**(1.0/3.0)
    
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

subroutine pbc3d (vec,L)

!Vector should be in the range 0 -> L in all dimensions
!Therefore, we need to apply the following changes if it's not in this range:
    implicit none
    real*8,dimension(3),intent(inout) :: vec
    real*8, intent(in) :: L
    integer :: i
    
    do i=1,3
        do while (vec(i).gt.L) 
            vec(i) = vec(i)-L
        end do
        
        do while (vec(i).lt.0) 
            vec(i) = vec(i)+L
        end do
    
    end do
    





end subroutine pbc3d

! calculate energy of lj fluid

real*8  function LJ_Energy(coords,L,npart,density)

    implicit none
    real*8, dimension(3,300),intent(in) :: coords
    real*8, intent(in) :: L, density
    integer, intent(in) :: npart
    integer :: partA,partB
    real*8,dimension(3) :: dr
    LJ_Energy = 0
    
    do partA = 1, npart-1
        do partB = partA+1,npart
        !calculate particle partical distance
        dr(1) = coords(1,partA)-coords(1,partB)
        dr(2) = coords(2,partA)-coords(2,partB)
        dr(3) = coords(3,partA)-coords(3,partB)
   
       call distpbc3d (dr,L) ! apply periodic bdry condition
        
      
        if((dot_product(dr,dr))**(0.5).lt.0.5*L) then
            LJ_Energy = LJ_Energy+(((1.0/dot_product(dr,dr))**(3))*(((1.0/dot_product(dr,dr))**(3))-1))
        end if
        
        end do
    
    
    
    end do
    
    LJ_Energy=LJ_Energy*4
    LJ_Energy=LJ_Energy+npart*(8*3.141592/(3.0))*density*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
    !tail correction


end function LJ_Energy

real*8  function LJ_EnergyChange(coords,rtrial,j,L,npart)

    implicit none
    real*8, dimension(3,300),intent(in) :: coords
    real*8, dimension(3),intent(in) :: rtrial
    real*8, intent(in) :: L
    integer, intent(in) :: npart,j
    integer :: otherPart
    real*8,dimension(3) :: drnew, drold
    LJ_EnergyChange = 0
    
    do otherPart = 1, npart
        
        if(j.ne.otherPart) then
        !calculate particle partical distance
        
        drnew(1) = coords(1,otherPart)-rtrial(1)
        drnew(2) = coords(2,otherPart)-rtrial(2)
        drnew(3) = coords(3,otherPart)-rtrial(3)
        
        drold(1) = coords(1,otherPart)-coords(1,j)
        drold(2) = coords(2,otherPart)-coords(2,j)
        drold(3) = coords(3,otherPart)-coords(3,j)
      
     call distpbc3d (drnew,L) ! apply periodic bdry condition
     call distpbc3d (drold,L) ! apply periodic bdry condition
        
            if((dot_product(drnew,drnew))**(0.5).lt.0.5*L) then
                LJ_EnergyChange = LJ_EnergyChange+&
                &(((1.0/dot_product(drnew,drnew))**(3))*(((1.0/dot_product(drnew,drnew)))**(3)-1))
                
            end if
            
            if((dot_product(drold,drold))**(0.5).lt.0.5*L) then 
                LJ_EnergyChange = LJ_EnergyChange &
                &-(((1.0/dot_product(drold,drold))**(3))*(((1.0/dot_product(drold,drold))**(3))-1))
            end if
        
        end if
    
    
    
    end do
    
    LJ_EnergyChange=LJ_EnergyChange*4
    



end function LJ_EnergyChange

subroutine distpbc3d (vec,L)

!Vector should be in the range 0 -> L in all dimensions
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
subroutine mcexc (L,V,npart,coords,zz,beta)

    implicit none
    real*8 :: densityo,arg,densityn
    real*8, intent(in) :: beta,V, L,zz
    integer, intent(inout) :: npart
    real*8, dimension(3,300), intent(inout) :: coords
    real*8, external :: LJ_EnergyChange1
    integer :: o
    real*8, dimension(3) :: xn
    real*8 :: enn, eno
    
   
    if (rand().lt.0.5) then
        densityo = npart/(L**3)
        densityn = (npart-1.0)/(L**3)
        o = int(npart*rand())+1 ! select a particle to be removed
        eno = LJ_EnergyChange1(coords,coords(:,o),o,L,npart) !energy particle o
        eno = eno - (npart-1.0)*(8*3.141592/(3.0))*densityn*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        eno = eno + (npart)*(8*3.141592/(3.0))*densityo*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        
        ! tail correction
        arg = npart*exp(beta*eno)/(zz*V)
        if (rand().lt.arg) then
            coords(:,o) = coords(:,npart)
            npart = npart - 1
        end if
    else
    densityo = npart/(L**3)
    densityn = (npart+1.0)/(L**3)
        xn(1) = rand()*L !new particle at a random position
        xn(2) = rand()*L
        xn(3) = rand()*L
        enn = LJ_EnergyChange1(coords,xn,200,L,npart)
        enn = enn + (npart+1)*(8*3.141592/(3.0))*densityn*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        enn = enn - (npart)*(8*3.141592/(3.0))*densityo*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        arg = zz*V*exp(-beta*enn)/(npart+1)
        if (rand().lt.arg) then
            coords(:,npart+1) = xn(:)
            npart = npart + 1
        end if
    end if
        
end subroutine mcexc

real*8  function LJ_EnergyChange1(coords,vec,j,L,npart)

    implicit none
    real*8, dimension(3,300),intent(in) :: coords
    real*8, dimension(3),intent(in) :: vec
    real*8, intent(in) :: L
    integer, intent(in) :: npart,j
    integer :: otherPart
    real*8,dimension(3) :: dr
   
    LJ_EnergyChange1 = 0
    
    do otherPart = 1, npart
        
        if(j.ne.otherPart) then
        !calculate particle partical distance
        
        dr(1) = coords(1,otherPart)-vec(1)
        dr(2) = coords(2,otherPart)-vec(2)
        dr(3) = coords(3,otherPart)-vec(3)
        
     call distpbc3d (dr,L) ! apply periodic bdry condition
        
            if((dot_product(dr,dr))**(0.5).lt.0.5*L) then
                LJ_EnergyChange1 = LJ_EnergyChange1+&
                &(((1.0/dot_product(dr,dr))**(3))*(((1.0/dot_product(dr,dr)))**(3)-1))
                
            end if
            
        
        end if
    
    
    
    end do
    
    LJ_EnergyChange1=LJ_EnergyChange1*4
   
   



end function LJ_EnergyChange1

