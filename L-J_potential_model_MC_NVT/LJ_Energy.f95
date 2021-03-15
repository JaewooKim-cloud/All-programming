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
        dr = coords(:,partA)-coords(:,partB)
       call pbc3d (dr,L) ! apply periodic bdry condition
        
        
        
        LJ_Energy = LJ_Energy+(((1.0/dot_product(dr,dr))**(3))*(((1.0/dot_product(dr,dr))**(3))-1))
        
        
        end do
    
    
    
    end do
    
    LJ_Energy=LJ_Energy*4
    



end function LJ_Energy
