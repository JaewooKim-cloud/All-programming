program mc_nvt

          ! we apply widom method
            
    implicit none
    integer,parameter :: seed = 810907
    integer, parameter :: npart = 100, nsteps = 20000, printFreq = 100
    real*8 :: density , temp, beta, maxDr ,rc,rc2,ecut
    real*8, dimension (3,npart) :: coords
    integer i,j
    real*8 :: L, energy,deltaE,wtest,muex
    real*8,dimension(3) :: rtrial
    real*8,external :: LJ_Energy
    real*8,external :: LJ_EnergyChange
    real*8,external :: LJ_EnergyChange1
    
    wtest = 0
    density = 0.7
    temp = 2.0
    beta = 1.0/temp
    maxDr = 0.4
    
    rc = 0.5*L
    rc2 = (rc)**2
    ecut = 4*(1.0/(rc**12)-1.0/(rc**6))
    call srand(seed)
    
   call initcubicgrid(npart,density,L,coords) ! set initial configuraion
   
   energy = LJ_Energy(coords,L, npart)
   
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
                energy = energy + deltaE
            end if
            if (i.gt.nsteps*0.5) then
          call  widom (npart,L,beta,coords,wtest)
            
        end if
        end do
        
   
    end do
    
    muex = -(1.0/beta)*log(wtest/(0.5*nsteps*npart))
    muex= muex + 2*(8*3.141592/(3.0))*density*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
       
   print *, energy
    print *, muex

end program mc_nvt


subroutine initcubicgrid(npart, density,L,coords)

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

subroutine pbc3d (vec,L)

!Vector should be in the range 0 -> L in all dimensions
!Therefore, we need to apply the following changes if it's not in this range:
    implicit none
    real*8,dimension(3),intent(inout) :: vec
    real*8, intent(in) :: L
    integer :: i
    
    do i=1,3
        if (vec(i).gt.L) then
            vec(i) = vec(i)-L
        end if
        
        if (vec(i).lt.0) then
            vec(i) = vec(i)+L
        end if
    
    end do
    





end subroutine pbc3d

! calculate energy of lj fluid

real*8  function LJ_Energy(coords,L,npart)

    implicit none
    real*8, dimension(3,npart),intent(in) :: coords
    real*8, intent(in) :: L
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
        
      
        
        LJ_Energy = LJ_Energy+(((1.0/dot_product(dr,dr))**(3))*(((1.0/dot_product(dr,dr))**(3))-1))
        
        
        end do
    
    
    
    end do
    
    LJ_Energy=LJ_Energy*4
    LJ_Energy=LJ_Energy+npart*(8*3.141592/(3.0))*npart*(1.0/(L**3))*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        



end function LJ_Energy

real*8  function LJ_EnergyChange(coords,rtrial,j,L,npart)

    implicit none
    real*8, dimension(3,npart),intent(in) :: coords
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
  &(((1.0/dot_product(drnew,drnew))**(3))*(((1.0/dot_product(drnew,drnew))**(3))-1)) &
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

subroutine widom (npart,L,beta,coords,wtest)!calculate excess chemical potential
        
        implicit none
        integer, intent(in) :: npart
        real*8 :: densityo,densityn,entest
        real*8, dimension(3) :: xtest
        real*8, intent(in) :: L,beta
        real*8,external :: LJ_EnergyChange1
        real*8,dimension(3,npart) ::coords
        real*8, intent(inout) :: wtest
        
        densityo = npart/(L**3)
        densityn = (npart+1.0)/(L**3)
        xtest(1) = L*rand()
        xtest(2) = L*rand()
        xtest(3) = L*rand()
        entest = LJ_EnergyChange1(coords,xtest,200,L,npart)
        ! tail correction
        !entest = entest + (npart+1)*(8*3.141592/(3.0))*densityn*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        !entest = entest - (npart)*(8*3.141592/(3.0))*densityo*((((1.0/(0.5*L))**9)/3.0-(1.0/(0.5*L))**3))
        
        wtest = wtest + exp(-beta*entest)

end subroutine widom
