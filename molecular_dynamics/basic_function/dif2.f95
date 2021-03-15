subroutine dif2(switch,nsamp,ntel,dtime,dt,ibl,npart,tel,delr2,vxsum,v,iblm,ibmax,n)

    implicit none
    integer, intent(in) :: switch, nsamp,npart,ibmax, n
    integer, intent(inout) :: ntel,iblm
    real*8, intent(inout) :: dtime
    real*8, intent(in) :: dt
    integer :: ib,i,j,ii
    integer, dimension (ibmax),intent(inout) :: ibl
    integer, dimension (ibmax,npart),intent(inout) :: tel
    real*8, dimension (ibmax,npart),intent(inout) :: delr2
    real*8, dimension (ibmax,n,npart),intent(inout) :: vxsum
    real*8 :: r2,delx,time
    real*8,parameter :: tdifmax = 2
    integer :: inm, inp, in
    real*8, dimension(3,npart), intent(in) :: v
    
    if (switch.eq.0) then
    
        ntel = 0
        dtime = dt*nsamp
        
        do ib = 1, ibmax
         ibl(ib) = 0 !length of current block
         do j=1,n
            tel (ib,j)=0
            delr2(ib,j)=0
            do i=1,npart
                vxsum(ib,j,i)=0 ! coarse-grained velocity particle i
            end do
         end do
        end do
    else if (switch.eq.2) then
       
       open(unit=1,file="diffusion2.txt")
       10 format(f10.4, 5x , f10.4,5x)
       
        do ib=1,min(ibmax,iblm)
            do j=2,min(ibl(ib),n)
                time=dtime*j*n**(ib-1)
                r2 = delr2(ib,j)*(dtime**2)/tel(ib,j)
                write(unit=1,fmt=10),time, r2
            end do
        end do
        
        close(unit=1)
    else if (switch.eq.1) then
        ntel = ntel + 1
        !determine current maximum number of blocks : iblm
        iblm = 1
        ii = ntel/n
        do while (ii.ne.0)
            iblm = iblm + 1
            ii= ii/n
            !test maximum time not longer than tdimax
            
            if (dtime*(n**iblm).gt.tdifmax) ii=0 !tdifmax : maximum time to determine diffusion
        end do
        iblm = min(iblm,ibmax) !normally, n=10 and ibmax = 20
        
        do ib = 1, iblm
            if (mod(ntel,n**(ib-1)).eq.0) then
                ibl(ib) = ibl(ib) + 1 !increase current block length
                inm = min(ibl(ib),n) !set maximum block length to n
                do i=1,npart
                    if (ib.eq.1) then
                        delx = v(1,i) !0th block : ordinary velocity
                    else
                        delx = vxsum(ib-1,1,i) !previous block velocity
                    end if
                    do in = 1,inm
                        if (inm.ne.n) then
                            inp = in
                        else
                            inp = in + 1
                        end if
                        if(in.lt.inm) then
                            vxsum(ib,in,i) = vxsum(ib,inp,i) + delx
                        else
                            vxsum(ib,in,i) = delx
                        end if
                    end do
                    do in=1,inm
                        tel(ib,in) = tel(ib,in) + 1 !counter number of updates
                        delr2(ib,in) = delr2(ib,in) + vxsum(ib,inm-in+1,i)**2
                    end do
                end do
            end if
        end do
    end if



end subroutine dif2
