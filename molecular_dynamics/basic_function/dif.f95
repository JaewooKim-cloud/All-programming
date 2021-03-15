subroutine dif (switch,nsamp,t0max,it0,tmax,ntel,t0,dtime,dt,ntime,vacf,r2t,time,x0,v0,x,v,npart)

    implicit none
    integer, intent(in) :: switch,nsamp,t0max, it0, tmax,npart
    integer, intent(inout) :: ntel,t0
    real*8, intent(inout) :: dtime,dt
    
    integer,dimension(tmax) ,intent(inout) ::ntime
    real*8,dimension(3,tmax) ,intent(inout) :: vacf
    real*8,dimension(tmax) ,intent(inout) :: r2t, time
    integer :: tt0,i,t,delt
    real*8,dimension(3,npart,t0max) ,intent(inout) :: x0,v0
    real*8,dimension(3,npart) ,intent(in) :: x,v
    integer, dimension(t0max) :: time0
    
    if (switch.eq.0) then ! initialization
        ntel = 0 !time counter
        dtime = dt*nsamp
        do i=1,tmax ! total number of time step
            ntime(i) = 0 ! number of samples for time i
            vacf(:,i) = 0
            r2t(i) = 0
        end do
    else if(switch.eq.1)then
        ntel = ntel + 1
        if (mod(ntel,it0).eq.0) then
            t0 = t0 + 1 !update number of t=o
            tt0 = mod(t0-1,t0max)+1
            time0(tt0) = ntel !store the time of t0
            do i=1, npart
                x0(:,i,tt0) = x(:,i)
                v0(:,i,tt0) = v(:,i)
            end do
        end if
        do t=1,min(t0,t0max)
            delt = ntel - time0(t) + 1
            if (delt.lt.tmax) then
                ntime(delt) = ntime(delt) + 1
                do i = 1,npart
                    vacf(:,delt) = vacf(:,delt) + v(:,i)*v0(:,i,t)
                    r2t(delt) = r2t(delt) + dot_product(x(:,i)-x0(:,i,t),x(:,i)-x0(:,i,t))
                end do
            end if
        end do
    else if (switch.eq.2) then
        do i = 1, tmax
            time(i) = dtime*(i+0.5)
            vacf(:,i) = vacf(:,i)/(npart*ntime(i))
            r2t(i) = r2t(i)/(npart*ntime(i))
        end do
    end if

end subroutine dif
