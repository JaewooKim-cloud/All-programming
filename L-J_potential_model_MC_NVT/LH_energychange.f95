

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
      
     
        
        
        LJ_EnergyChange = LJ_EnergyChange+&
  &(((1.0/dot_product(drnew,drnew))**(3))*(((1.0/dot_product(drnew,drnew))**(3))-1)) &
  &-(((1.0/dot_product(drold,drold))**(3))*(((1.0/dot_product(drold,drold))**(3))-1))
        
        
        end if
    
    
    
    end do
    
    LJ_EnergyChange=LJ_EnergyChange*4
    



end function LJ_EnergyChange
