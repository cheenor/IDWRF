subroutine getcomdbz(dbz,nz,ny,nx,comdbz)
    implicit none
    integer,intent(in) :: nz,ny,nx
    real, dimension(nz,ny,nx), intent(in):: dbz
    real, dimension(ny,nx), intent(out):: comdbz
    !
    integer ix,iy,iz
    real temp
    do ix=1,nx
        do iy=1,ny
            temp=-9999
            comdbz(iy,ix)=0
            do iz=1,nz
                if(dbz(iz,iy,ix)>temp)then
                    temp=dbz(iz,iy,ix)
                endif
            enddo
            comdbz(iy,ix)=temp
        enddo       
    enddo
    return
end subroutine getcomdbz
subroutine get_mid_xz(ain,ny1,ny2,nz,ny,nx,aout)
    implicit none
    integer, intent(in):: ny1,ny2
    integer, intent(in):: nz,ny,nx
    real, dimension(nz,ny,nx),intent(in) :: ain
    real, dimension(nz,nx),intent(out) :: aout
    !
    integer iz,iy,ix
    real temp,a
    a=1.0*(ny2-ny1+1)
    do iz=1,nz
        do ix=1,nx
            aout(iz,ix)=0.0
            do iy=ny1,ny2
                aout(iz,ix)=aout(iz,ix)+ain(iz,iy,ix)/a
            enddo
        enddo
    enddo
    return 
end subroutine get_mid_xz
subroutine domain_prf(ain,nz,ny,nx,aout)
    implicit none
    integer, intent(in):: nz,ny,nx
    real, dimension(nz,ny,nx),intent(in) :: ain
    real, dimension(nz),intent(out) :: aout
    !
    integer iz,iy,ix
    do iz=1,nz
        aout(iz)=0.0
        do ix=1,nx
            do iy=1,ny
                aout(iz)=aout(iz)+ain(iz,iy,ix)
            enddo
        enddo
        aout(iz)=aout(iz)/(nx*ny*1.0)
    enddo
    return 
end subroutine domain_prf
subroutine domain_2d_prf_lwc(qc,qr,nz,ny,nx,aout)
    implicit none
    integer, intent(in):: nz,ny,nx
    real, dimension(nz,ny,nx),intent(in) :: qc,qr
    real, dimension(nz),intent(out) :: aout
    !
    real lwc_bin(10)
    integer iz,iy,ix
    do iz=1,nz
        aout(iz)=0.0
        do ix=1,nx
            do iy=1,ny
                aout(iz)=aout(iz)+qc(iz,iy,ix)+qr(iz,iy,ix)
            enddo
        enddo
        aout(iz)=aout(iz)/(nx*ny*1.0)
    enddo
    return 
end subroutine domain_2d_prf_lwc
subroutine get_cell_chars(qc,qr,nz,ny,nx,aout)
    implicit none
    integer, intent(in):: nz,ny,nx
    real, dimension(nz,ny,nx),intent(in) :: qc,qr
    real, dimension(nx),intent(out) :: aout
    !
    integer iz,iy,ix,i
    real tcld(nz)
    integer cb(ny,nx,nz),ct(ny,nx,nz)
    integer ncell,k1,k2
    real aa,cm(ny,nx,nz),dmx,cmx
    integer dmx_bt(2),cmx_bt(2),dpt
    integer nct(ny,nx)
    cb=0
    ct=0
    nct=0
    do ix=1,nx
        do iy=1,iy
            do iz=1,nz
                tcld(iz)=qc(iz,iy,ix)*1000.+qr(iz,iy,ix)*1000. ! convert to g kg             
            enddo
            !# detaching the cloud cells
            ncell=0
            if (tcld(1)>0.01)then
                k1=1
                k2=2
                aa=tcld(1)
            else
                aa=0
            endif
            do iz=2,nz
                IF (tcld(iz-1).LE.0.01.AND.tcld(iz).GT.0.01) THEN
                    K1 = iz
                    K2 = iz
                    AA = MAX(AA,tcld(iz))
                ELSEIF (tcld(iz-1).GT.0.01.AND.tcld(iz).GT.0.01) THEN
                    K2 = iz
                    AA = MAX(AA,tcld(iz))
                ELSEIF (tcld(iz-1).GT.0.01.AND.tcld(iz).LE.0.01) THEN
                    ncell = ncell + 1
                    cb(iy,ix,ncell) = K1
                    ct(iy,ix,ncell) = K2
                    CM(iy,ix,ncell) = AA
                    AA = 0.
                ENDIF
            enddo
            nct(iy,ix)=ncell
        enddo
    enddo
    dmx=0.0
    cmx=0.0
    do ix =1,nx
        do iy=1,ny
            ncell=nct(iy,ix)
            do i=1,ncell
                aa=cm(iy,ix,i)
                dpt=ct(iy,ix,i)-cb(iy,ix,i) 
                if (aa>cmx)then
                    cmx_bt(1)=cb(iy,ix,i)
                    cmx_bt(2)=ct(iy,ix,i)
                    cmx=aa
                endif
                if(dpt>dmx)then
                    dmx_bt(1)=cb(iy,ix,i)
                    dmx_bt(2)=ct(iy,ix,i)
                    dmx=dpt
                endif
            enddo
        enddo
    enddo
    aout(1)=cmx_bt(1)
    aout(2)=cmx_bt(2)
    aout(3)=cmx
    aout(4)=dmx_bt(1)
    aout(5)=dmx_bt(2)
    aout(6)=dmx
    return 
end subroutine get_cell_chars


