subroutine calmean(dbz,rainnc,rainc,ny1,ny2,ny,nx,rdbz,rnnc,rnc)
    implicit none
    integer,intent(in) :: ny,nx
    integer,intent(in) :: ny1,ny2
    real, dimension(ny,nx), intent(in):: dbz,rainnc,rainc
    real, dimension(nx), intent(out):: rdbz,rnnc,rnc
    !
    integer ix,iy
    real asum1,asum2,asum3,y,cont
    y=1.0/(ny*1.0)
    do ix=1,nx
        asum1=0.0
        asum2=0.0
        asum3=0.0
        cont=0.0
        do iy=ny1,ny2
            asum1=asum1+rainnc(iy,ix)
            asum2=asum2+rainc(iy,ix)
            if(dbz(iy,ix)>-30)then
            asum3=asum3+dbz(iy,ix)
            cont=cont+1.0
            endif
        enddo
        rnnc(ix)=asum1/(ny*1.0)
        rnc(ix)=asum2/(ny*1.0)
        if(cont>0)then
        rdbz(ix)=asum3/cont
        else
        rdbz(ix)=0.0
        endif
    enddo
    return
end subroutine calmean
subroutine calmeanv3(comdbz,rainnc,lats,lons,ny,nx,rdbz,rnnc)
     implicit none
    integer,intent(in) :: ny,nx
    real, dimension(ny,nx), intent(in):: comdbz,rainnc,lats,lons
    real, dimension(nx), intent(out):: rdbz,rnnc
    !
    integer ix,iy
    real asum1,asum2,asum3,y,cont,cont1
    y=1.0/(ny*1.0)
    rdbz=0.0
    rnnc=0.0
    do ix=1,nx
        asum1=0.0
        asum2=0.0
        asum3=0.0
        cont=0.0
        cont1=0.0
        do iy=1,ny  !ny1,ny2
            if (lons(iy,ix)>119 .and. lons(iy,ix)<123) then
                if (lats(iy,ix)>29 .and. lats(iy,ix)<33) then
                    asum1=asum1+rainnc(iy,ix)
                    cont1=cont1+1.0
                    if(comdbz(iy,ix)>0)then
                    asum3=asum3+comdbz(iy,ix)
                    cont=cont+1.0
                    endif
                endif
            endif
        enddo
        if (cont1>0)then
        rnnc(ix)=asum1/cont1
        endif
        if(cont>0)then
        rdbz(ix)=asum3/cont
        else
        rdbz(ix)=-99
        endif
    enddo
    return
end subroutine calmeanv3   
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
subroutine getbndlat(lats,blat,elat,ny,ny1,ny2)
    implicit none
    integer,intent(in) :: ny
    real, dimension(ny), intent(in):: lats
    real,  intent(in):: blat,elat
    integer, intent(out):: ny1,ny2
    integer ix,iy
    real aa0,aa1
    ny1=0
    ny2=0
    do iy=2,ny
        aa0=lats(iy-1)
        aa1=lats(iy)
        if(aa0<=blat .and. aa1>blat)then
            ny1=iy
        endif
        if(aa0<=elat .and. aa1>elat)then
            ny2=iy
        endif
    enddo       
    return
end subroutine getbndlat
subroutine regrid(lons,lats,rainnc,tlons,tlats,nx,ny,nyy,nxx,arainnc)
    implicit none
    integer,intent(in) :: ny,nx
    integer,intent(in) :: nxx,nyy
    real, dimension(ny), intent(in):: tlats
    real, dimension(nx), intent(in):: tlons
    real, dimension(nyy,nxx), intent(in):: lons,lats,rainnc
    real, dimension(ny,nx), intent(out):: arainnc
    !
    integer i,j,ii,jj,is,ie,js,je
    integer ix,iy,cont,idata
    integer ixx,iyy,m,n
    real*8 a,b,aa,bb,ds
    real*8 r,pi,dlh,dlv,dh,dv
    real*8   cri,td,tmda
    real   defv,minv,maxv,ser,dg,dsc
    integer itar,iser,iii,jjj
    m=nyy
    n=nxx
    !
    minv=0.
    maxv=99999.
    pi=3.141592657
    r=6371.393
    dlv=2*pi*r/360.0
    cri=15. ! 5km
    iser=floor(5*3./3.)+5  ! 
    is=1;ie=m
    js=1;je=n
    idata=0
    defv=0
    arainnc=defv
    ii=0;jj=0
    iii=0;jjj=0
    do ix=1,nx
        a=tlons(ix)
!       ii=0
        ii=0;jj=0
        iii=0;jjj=0
        do iy=1,ny
            b=tlats(iy)
            dlh=2*pi*r*cos(b*pi/180.)/360.
            td=0
            tmda=0
            itar=0
            is=1;ie=m
            js=1;je=n
            if (ii>0 .and. jj>0)then
                if (iy==1)then
                    if(ix==1)then
                        iii=ii;jjj=jj
                    elseif(ix>1)then
                        ii=iii; jj=jjj
                        iii=ii ; jjj=jj
                    endif
                endif
                is=ii-iser; ie=ii+iser
                if(is<1)is=1
                if(ie>m)ie=m
                js=jj-iser; je=jj+iser 
                if(js<1)js=1
                if(je>n)je=n 
!               print*,is,ie,js,je
!               pause
            endif
            dsc=0.
!           print*,is,ie,js,je
            ii=0;jj=0
            arainnc(ix,iy)=0
            do i=1,nyy !is,ie
                do j=1,nxx !js,je
                    aa=lons(i,j) ; bb=lats(i,j)
                    dh=abs(aa-a)*dlv
                    if (dh<cri)then
                        dv=abs(bb-b)*dlh
                        ds=dh*dh+dv*dv
                        if (ds<cri*cri)then
                            ii=i  ;  jj=j
 !                          print*,ii,jj
                            if (dsc==0.)dsc=ds
                            if(rainnc(i,j)>= minv  &
     &                            .and. rainnc(i,j)<= maxv)then
                                if(dsc>=ds)then
                                    ii=i;jj=j
                                    dsc=ds
                                endif  
                                if (ds>0.)then
                                    td=td+1./ds                                               
                                    tmda=tmda+rainnc(i,j)/ds
                                elseif(ds==0)then
                                    arainnc(ix,iy)=rainnc(i,j)
                                    itar=1
                                    exit
                                endif
                            endif
                        endif
                    endif
                enddo
                if(itar==1)then
                    exit
                endif
            enddo 
            if (td/=0. .and. itar==0)then 
                arainnc(ix,iy)=tmda/td
                if(arainnc(ix,iy)>0)then
                    print*,arainnc(ix,iy)
                endif
            endif                
        enddo      
    enddo
    return
end subroutine regrid
subroutine calmeanv2(comdbz,rainnc,ny1,ny2,ny,nx,rdbz,rnnc)
    implicit none
    integer,intent(in) :: ny,nx
    integer,intent(in) :: ny1,ny2
    real, dimension(ny,nx), intent(in):: comdbz,rainnc
    real, dimension(nx), intent(out):: rdbz,rnnc
    !
    integer ix,iy
    real asum1,asum2,asum3,y,cont
    y=1.0/((ny2-ny1+1.0)*1.0)
    do ix=1,nx
        asum1=0.0
        asum2=0.0
        asum3=0.0
        cont=0.0
        do iy=ny1,ny2
            asum1=asum1+rainnc(iy,ix)
            if(comdbz(iy,ix)>-30)then
                asum3=asum3+comdbz(iy,ix)
                cont=cont+1.0
            endif
        enddo
        rnnc(ix)=asum1*y
        if(cont>0)then
            rdbz(ix)=asum3/cont
        else
            rdbz(ix)=0.0
        endif
    enddo
    return
end subroutine calmeanv2