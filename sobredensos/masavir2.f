c-------gfortran -mcmodel=medium -o masavir2 masavir2.f  
c       (--> ya esta en el alias gfortran)
c       -mcmodel=medium --> amplia el uso de ram
c       -------------------------------------------- 

!gfortran -mcmodel=large -o masavir2 masavir2.f
!input1: cosmologia 1, 2 o 3
!input2: redshift 0-27
!./masavir2 input1 input2     
c--------------------
c La magnitud absoluata NO viene con el -5log(h)   
c--------------------         
        program masavir2
        implicit none
        integer i,j,k,h,indr,igrumax,ngal,ngaln,ii
        integer ngfin
        integer snap,ncosmo  
        integer nmax,imax,gmax,nngrid,nbin
        parameter (imax=15000,gmax=9000000,nngrid=128,nbin=15)
        parameter (nmax=100000000) !100 millones
        real den_num_med
        real G,c
        real h0
c ------galaxias
        integer first(gmax),ll(nmax),mark(nmax),nrep
        integer lirst(nngrid,nngrid,nngrid),lg(nmax),nvec
        integer lirstn(nngrid,nngrid,nngrid),lgn(nmax)
        real dmin,xvec,yvec,zvec
        real xbin(nbin),dens_cmp(nbin),dens_new(nbin),dens_mie(nbin),dens_gru(nbin)
        real den_aisl,den_a,cascaron,pi,den_cmp

        integer igru(nmax),inum(nmax),nmie(gmax),iigru(imax),iinum(imax)
        integer*8 gID(nmax),ggID(imax)
        integer it,igrut,inumt
        integer igrut_df,inumt_df,igru_df(nmax),inum_df(nmax),iigru_df(imax),iinum_df(imax)
        integer*8 gIDt
        real umagt,gmagt,rmagt,imagt,zmagt
        real kmagt,jmagt !magnitudes del Hen_15
        real stmasst,bmasst
        real*8 dist
        real umag(nmax),gmag(nmax),rmag(nmax),imag(nmax),zmag(nmax)
        real kmag(nmax),jmag(nmax)
        real stmass(nmax),bmass(nmax),munit
        real x,y,z,difx,dify,difz
        real vx,vy,vz,vp
        real xt(nmax),yt(nmax),zt(nmax)
        real vxt(nmax),vyt(nmax),vzt(nmax)
        real zct(nmax),zst(nmax)
        real rr_new(3,nmax),rmag_new(nmax),stm_new(nmax)
        integer igru_new(nmax)
        real fac,rl,rsize,rr(3,nmax)
        integer ngrid
        integer ngal_s
c-------galaxias en grupos
        real xx(imax),yy(imax),zz(imax),lr(imax)
        real*8 rrmag(imax)
        real ggmag(imax),sstmass(imax),bbmass(imax)
        real uumag(imax),iimag(imax),zzmag(imax)
        real kkmag(imax),jjmag(imax)
        real xxcm,yycm,zzcm
        real xcm,ycm,zcm
        real xx2(imax),yy2(imax),zz2(imax)
        real xcm1
c-------galaxias en grupos
        real vvx(imax),vvy(imax),vvz(imax)
        real*8 vv(imax),vvels(imax)
        real vvxcm,vvycm,vvzcm
        real vxcm,vycm,vzcm
        real zzz,zc(imax),zs(imax)
        real v(nmax),vels(nmax),vels_max
        real rvir_3d,sig_3d,rvir_3d_old
        real rvir_bi,rij(nmax)
        integer n_rij
        real mass_bi,lbox,lb2
        real numx,numy,numz,denx,deny,denz,xl,yl,zl
        real lsun, msun
        real rabs1, rabs_n, deltamag
        real zc_min, zc_max, delta_zc
        real zs_min, zs_max, delta_zs
        integer indx(imax),indc(imax),inds(imax)
        integer n3,nd
        real xij,yij,zij
        real*8 dist_ij(nmax)
        real dist_media
        real dist_mediana,sig_med
        real dmax,dmax_x,suma_stm
        integer indsat(imax),im, indv(imax)
        real mst_1,mst_2,mst_3,coef_mst,mst_1_bri

c---------------------------------------------------------------------
        external  iargc
        integer iargc
        character*80 string
        character*80 filein,filecosmo
        character*80 fileout1,fileout2,fileout3,fileout4

c       Tabla redshift
        integer ntab
        real ztab(70000)
        real rtab(70000)
        common /tabla/ ztab,rtab,ntab
c   
        integer ipcall,iplast,iphead,jc
        common /kpercent/ ipcall,iplast,iphead,jc

        if (iargc().ne.2) then
        print *,'usage:  ncosmo snap'
        print *,'Wmap_1/guo11 1 '
        print *,'Wmap_7/guo13 2 '
        print *,'Planck/hen15 3 '        
        print *,'wmap1/guo_II 4 '        
        print *,'wmap1/guo_II_16 5 '        
        print *,'Planck/hen20 7 '        
        print *,'Planck/ayr21 8 '        
        print *,'snap: z en el que estas parada (ejemplo: 5)'
        stop
        endif
        

        call getarg(1,string)
        read(string,*)ncosmo 
        call getarg(2,string)
        read(string,*)snap             

        c=299792.458 !Km/s
        G=6.6726*10.**(-11.)  !m**3/Kg s**2
        munit=(3.0857*10.**(22.))/(1.99*10.**(30.))
        lsun=3.846*(10.**33.) !erg/s
        msun=4.65 !r-band http://mips.as.arizona.edu/~cnaw/sun.html 
        pi=4.*atan(1.)

c------------ FILENAMES -------------------------------   
        if(ncosmo.eq.1)then
          lbox=500.
          h0=0.73
          filecosmo='Mill_I/red2dis_lcdm.dat-mill_wmap1'
          if(snap.lt.10)then
           write(filein,'("guo_11/z",i1,"/tablaglx.dat")')snap
           write(fileout1,'("guo_11/z",i1,"/galx_gru.dat")')snap
           write(fileout2,'("guo_11/z",i1,"/tablagru.dat")')snap
           write(fileout3,'("guo_11/z",i1,"/perfiles.dat")')snap
           write(fileout4,'("guo_11/z",i1,"/file_dim.dat")')snap
        else
           write(filein,'("guo_11/z",i2,"/tablaglx.dat")')snap
           write(fileout1,'("guo_11/z",i2,"/galx_gru.dat")')snap
           write(fileout2,'("guo_11/z",i2,"/tablagru.dat")')snap
           write(fileout3,'("guo_11/z",i2,"/perfiles.dat")')snap
           write(fileout4,'("guo_11/z",i2,"/file_dim.dat")')snap
         endif
        endif
        !----------------
        if(ncosmo.eq.2)then
          lbox=500.
          h0=0.704
          filecosmo='Mill_I/red2dis_lcdm.dat-mill_wmap7'
          if(snap.lt.10)then
           write(filein,'("guo_13/z",i1,"/tablaglx.dat")')snap
           write(fileout1,'("guo_13/z",i1,"/galx_gru.dat")')snap
           write(fileout2,'("guo_13/z",i1,"/tablagru.dat")')snap
           write(fileout3,'("guo_13/z",i1,"/perfiles.dat")')snap
           write(fileout4,'("guo_13/z",i1,"/file_dim.dat")')snap
        else
           write(filein,'("guo_13/z",i2,"/tablaglx.dat")')snap
           write(fileout1,'("guo_13/z",i2,"/galx_gru.dat")')snap
           write(fileout2,'("guo_13/z",i2,"/tablagru.dat")')snap
           write(fileout3,'("guo_13/z",i2,"/perfiles.dat")')snap
           write(fileout4,'("guo_13/z",i2,"/file_dim.dat")')snap
         endif
        endif
        !----------------
        if(ncosmo.eq.3)then
          lbox=480.279
          h0=0.673
          filecosmo='Mill_I/red2dis_lcdm.dat-mill_planck'
          if(snap.lt.10)then
           write(filein,'("hen_15/z",i1,"/tablaglx.dat")')snap
           write(fileout1,'("hen_15/z",i1,"/galx_gru.dat")')snap
           write(fileout2,'("hen_15/z",i1,"/tablagru.dat")')snap
           write(fileout3,'("hen_15/z",i1,"/perfiles.dat")')snap
           write(fileout4,'("hen_15/z",i1,"/file_dim.dat")')snap
        else
           write(filein,'("hen_15/z",i2,"/tablaglx.dat")')snap
           write(fileout1,'("hen_15/z",i2,"/galx_gru.dat")')snap
           write(fileout2,'("hen_15/z",i2,"/tablagru.dat")')snap
           write(fileout3,'("hen_15/z",i2,"/perfiles.dat")')snap
           write(fileout4,'("hen_15/z",i2,"/file_dim.dat")')snap
         endif
        endif
        !----------------
       if(ncosmo.eq.4)then
          lbox=100.
          h0=0.73
          filecosmo='Mill_II/red2dis_lcdm.dat-mill_wmap1'
          if(snap.lt.10)then
           write(filein,'("guo_II/z",i1,"/tablaglx.dat")')snap
           write(fileout1,'("guo_II/z",i1,"/galx_gru.dat")')snap
           write(fileout2,'("guo_II/z",i1,"/tablagru.dat")')snap
           write(fileout3,'("guo_II/z",i1,"/perfiles.dat")')snap
           write(fileout4,'("guo_II/z",i1,"/file_dim.dat")')snap
        else
           write(filein,'("guo_II/z",i2,"/tablaglx.dat")')snap
           write(fileout1,'("guo_II/z",i2,"/galx_gru.dat")')snap
           write(fileout2,'("guo_II/z",i2,"/tablagru.dat")')snap
           write(fileout3,'("guo_II/z",i2,"/perfiles.dat")')snap
           write(fileout4,'("guo_II/z",i2,"/file_dim.dat")')snap
         endif
        endif
        !----------------
       if(ncosmo.eq.5)then
          lbox=100.
          h0=0.73
          filecosmo='Mill_II/red2dis_lcdm.dat-mill_wmap1'
          if(snap.lt.10)then
           write(filein,'("guo_II/z",i1,"_16/tablaglx.dat")')snap
           write(fileout1,'("guo_II/z",i1,"_16/galx_gru.dat")')snap
           write(fileout2,'("guo_II/z",i1,"_16/tablagru.dat")')snap
           write(fileout3,'("guo_II/z",i1,"_16/perfiles.dat")')snap
           write(fileout4,'("guo_II/z",i1,"_16/file_dim.dat")')snap
        else
           write(filein,'("guo_II/z",i2,"_16/tablaglx.dat")')snap
           write(fileout1,'("guo_II/z",i2,"_16/galx_gru.dat")')snap
           write(fileout2,'("guo_II/z",i2,"_16/tablagru.dat")')snap
           write(fileout3,'("guo_II/z",i2,"_16/perfiles.dat")')snap
           write(fileout4,'("guo_II/z",i2,"_16/file_dim.dat")')snap
         endif
        endif
        !----------------
        if(ncosmo.eq.7)then
          lbox=480.279
          h0=0.673
          filecosmo='Mill_I/red2dis_lcdm.dat-mill_planck'
          if(snap.lt.10)then
           write(filein,'("hen_20/z",i1,"/tablaglx.dat")')snap
           write(fileout1,'("hen_20/z",i1,"/galx_gru.dat")')snap
           write(fileout2,'("hen_20/z",i1,"/tablagru.dat")')snap
           write(fileout3,'("hen_20/z",i1,"/perfiles.dat")')snap
           write(fileout4,'("hen_20/z",i1,"/file_dim.dat")')snap
        else
           write(filein,'("hen_20/z",i2,"/tablaglx.dat")')snap
           write(fileout1,'("hen_20/z",i2,"/galx_gru.dat")')snap
           write(fileout2,'("hen_20/z",i2,"/tablagru.dat")')snap
           write(fileout3,'("hen_20/z",i2,"/perfiles.dat")')snap
           write(fileout4,'("hen_20/z",i2,"/file_dim.dat")')snap
         endif
        endif
        !----------------
        if(ncosmo.eq.8)then
          lbox=480.279
          h0=0.673
          filecosmo='Mill_I/red2dis_lcdm.dat-mill_planck'
          if(snap.lt.10)then
           write(filein,'("ayr_21/z",i1,"/tablaglx.dat")')snap
           write(fileout1,'("ayr_21/z",i1,"/galx_gru.dat")')snap
           write(fileout2,'("ayr_21/z",i1,"/tablagru.dat")')snap
           write(fileout3,'("ayr_21/z",i1,"/perfiles.dat")')snap
           write(fileout4,'("ayr_21/z",i1,"/file_dim.dat")')snap
        else
           write(filein,'("ayr_21/z",i2,"/tablaglx.dat")')snap
           write(fileout1,'("ayr_21hen_20/z",i2,"/galx_gru.dat")')snap
           write(fileout2,'("ayr_21/z",i2,"/tablagru.dat")')snap
           write(fileout3,'("ayr_21/z",i2,"/perfiles.dat")')snap
           write(fileout4,'("ayr_21/z",i2,"/file_dim.dat")')snap
         endif
        endif
        !----------------

c-------leo la tabla distancia vs redshift  ---------------------
        open(89,file=filecosmo,status='old')
        do i=1,nmax
                read(89,*,end=10)ztab(i),rtab(i)
        end do
 10     ntab=i-1
        write(*,*)'tabla de dist',ntab
        close(89)

        open(88,file=filein,status='old')
        open(90,file=fileout1,status='unknown')
        open(99,file=fileout2,status='unknown')
        open(33,file=fileout3,status='unknown')
        open(34,file=fileout4,status='unknown')
        
! Lectura de la tabla1
        !Inicializo todo
        igrumax=0
        ngal_s=0
        nmie=0
        xt=0.
        yt=0.
        zt=0.
        vp=0.
        it=0
        first=0
        ll=0
        suma_stm=0.


        lb2=lbox/2.
        rsize=10. !Mpc/h
        ngrid=int(lbox/rsize)
        rl=float(ngrid)
        fac=rl/lbox
        im=0
        do i=1,nmax
                call percent(i,13000000,'lectura')
                
                if(ncosmo.eq.3.or.ncosmo.eq.7.or.ncosmo.eq.8)then          
                          read(88,24,end=11)gIDt,x,y,z,dist,vx,vy,vz,
     &                    umagt,gmagt,rmagt,imagt,zmagt,kmagt,jmagt,
     &                    stmasst,bmasst,igrut,inumt,igrut_df,inumt_df
                else
                          read(88,23,end=11)gIDt,x,y,z,dist,vx,vy,vz,
     &                    umagt,gmagt,rmagt,imagt,zmagt,
     &                    stmasst,bmasst,igrut,inumt,igrut_df,inumt_df
                endif
              
c               if(rmagt.le.-16.)then   
                 !perfil de la GII c/todas las glxs (ncosmo=5)
                 im=im+1       
                 suma_stm=suma_stm+stmasst
                 !guardo pos,rmag,igru para glx en gru + campo para perfil
                 igru_new(im)=igrut
                 rr_new(1,im)=x*fac
                 rr_new(2,im)=y*fac
                 rr_new(3,im)=z*fac
                 rmag_new(im)=rmagt
                 stm_new(im)=stmasst
c                endif

                 !Contador solo glx en grup sobredensos
                 if(igrut.eq.0)ngal_s=ngal_s+1

                 if(inumt.lt.3)igrut=0 !gru glx menos 3 mie son campo
                 !Contador solo glx en grup sobredensos nmi>=3
                 if(igrut.eq.0)goto 21
                 nmie(igrut)=nmie(igrut)+1
                 if(igrut.GT.igrumax)igrumax=igrut
                       
                 it=it+1        
                 xt(it)=x
                 yt(it)=y
                 zt(it)=z
                 vxt(it)=vx
                 vyt(it)=vy
                 vzt(it)=vz
                 igru(it)=igrut
                 inum(it)=inumt
                 gID(it)=gIDt
                 !Magnitud absoluta
                 umag(it)=umagt-5.*log10(h0)
                 gmag(it)=gmagt-5.*log10(h0)
                 rmag(it)=rmagt-5.*log10(h0)
                 imag(it)=imagt-5.*log10(h0)
                 zmag(it)=zmagt-5.*log10(h0)
                 if(ncosmo.eq.3.or.ncosmo.eq.7.or.ncosmo.eq.8)then
                     kmag(it)=kmagt-5.*log10(h0) !Hen_15 y 20
                     jmag(it)=jmagt-5.*log10(h0) !Hen_15 y 20
                endif
                 stmass(it)=stmasst
                 bmass(it)=bmasst
                 igru_df(it)=igrut_df
                 inum_df(it)=inumt_df

                 !Lectura de la tabla de redshift
                 call locate(rtab,ntab,dist,indr)
                 zzz=(ztab(indr+1)-ztab(indr))/(rtab(indr+1)-
     &               rtab(indr))
                 zct(it)=zzz*(dist-rtab(indr))+ztab(indr)  

                 !velocidades pec=r.v/|r|
                 vp=(vx*x+vy*y+vz*z)/dist
                 zst(it)=(1.+zct(it))*(1.+vp/c)-1. 

                 v(it)=zst(it)*c

                 !velocidad gal
                 vels(it)=sqrt(vx**2+vy**2+vz**2)

                 first(igrut)=it
                 
                 rr(1,it)=x*fac
                 rr(2,it)=y*fac
                 rr(3,it)=z*fac
                 if(rr(1,it).gt.ngrid.or.rr(1,it).lt.0)write(*,*)'cagax',i
                 if(rr(2,it).gt.ngrid.or.rr(2,it).lt.0)write(*,*)'cagay',i
                 if(rr(3,it).gt.ngrid.or.rr(3,it).lt.0)write(*,*)'cagaz',i
c
 21     end do
 11     continue
        write(*,*)'numero de grupos sobredensos',igrumax
        ngal=it
        ngaln=im
        den_num_med=real(ngaln)/lbox**3
        
        write(*,*)'numero de galaxias en el box',ngaln
        write(*,*)'numero de galaxias en grupos sobredensos',ngal_s
        write(*,*)'numero de galaxias en grupos sobredensos nmi>=3',ngal
        write(*,*)'densidad (masa_st)',suma_stm/lbox**3
        write(*,*)'densidad numerica ',den_num_med
        write(34,*)igrumax,ngaln,den_num_med
        close(88)
c---------------------        
        do i=1,ngal
            ll(first(igru(i)))=i
            first(igru(i))=i 
         end do
c---------------------        
        lirst=0
        lg=0
        call grid(rr,ngrid,ngal,lirst,lg)
        lirstn=0
        lgn=0
        call grid(rr_new,ngrid,ngaln,lirstn,lgn)

c-----------para la percent--------------------------------
        ipcall=0
        iplast=0
        iphead=0
        jc=0
c---------------------        
        ! CALCULO PARA CADA GRUPO
        xx=0.
        yy=0.
        zz=0.
        vv=0.

        nrep=0
        mark=0
        ngfin=0
        do j=1,igrumax
               call percent(j,igrumax,'calculando')
               !write(*,*)j
c             restriccion en el numero de miembros para z altos
c             if(snap.eq.0)then
c               if (nmie(j).LT.10) goto 12
c             else
               if (nmie(j).LT.3) goto 12
c             end if
               ngfin=ngfin+1
               ii=0
               vvxcm=0.
               vvycm=0.
               vvzcm=0.
               i=first(j)
 87            i=ll(i)               
                        ii=ii+1
                        xx(ii)=xt(i)
                        yy(ii)=yt(i)
                        zz(ii)=zt(i)
                        rrmag(ii)=rmag(i)
                        ggmag(ii)=gmag(i)
                        sstmass(ii)=stmass(i)
                        bbmass(ii)=bmass(i)     
                        uumag(ii)=umag(i)
                        iimag(ii)=imag(i)
                        zzmag(ii)=zmag(i)
                        if(ncosmo.eq.3.or.ncosmo.eq.7.or.ncosmo.eq.8)then
                            kkmag(ii)=kmag(i) !Hen_15 y 20
                            jjmag(ii)=jmag(i) !Hen_15 y 20
                        endif
                        ggID(ii)=gID(i)
                        iigru(ii)=igru(i)
                        iinum(ii)=inum(i)
                        iigru_df(ii)=igru_df(i)
                        iinum_df(ii)=inum_df(i)
                        vvx(ii)=vxt(i)
                        vvy(ii)=vyt(i)
                        vvz(ii)=vzt(i)
                        zc(ii)=zct(i)
                        zs(ii)=zst(i)
                        vv(ii)=v(i)
                        vvels(ii)=vels(i)

                        !LUMINOSIDADES
                        lr(ii)=10.**(0.4*(msun-rmag(i)))

                        !VELOCIDAD PROMEDIO DE CADA GRUPO
                        vvxcm=vvxcm+vvx(ii)
                        vvycm=vvycm+vvy(ii)
                        vvzcm=vvzcm+vvz(ii)                        
               if(i.eq.first(j))goto 83
               goto 87
 83            continue

               if(ii.ne.nmie(j))write(*,*)'cuidado!!!'
               if(ii.ne.iinum(1))write(*,*)'cuidado2!!!',ii,
     1          iinum(1),j,nmie(j)
                vxcm=vvxcm/real(nmie(j))
                vycm=vvycm/real(nmie(j))
                vzcm=vvzcm/real(nmie(j))


                !----Calculo del radio virial biweight--------
                call dist_rij(nmie(j),lbox,xx,yy,zz,n_rij,rij)
                call rvir(1,n_rij,rij,rvir_bi)

                !-- clasico ------------------------------------
                !call rvir3d_old(nmie(j),lbox,xx,yy,zz,rvir_3d_old)
                !-- modificado ------------------------------------
                !call rvir3d(nmie(j),lbox,xx,yy,zz,sstmass,rvir_3d
                !if(rvir_3d.eq.0.)nnn=nnn+1
                !if(rvir_3d.eq.-99.)nnnn=nnnn+1
                !if(rvir_3d.eq.0.or.rvir_3d.eq.-99.)goto 12
                !------------------------------------

                call sig3d(nmie(j),vvx,vvy,vvz,vxcm,vycm,vzcm,sig_3d)
                mass_bi=(1000.*sig_3d)**2*rvir_bi/G*munit

                call indexx(ii,rrmag,indx)
                rabs1=rrmag(indx(1))   !brightest galaxy
                !rabs_n=rrmag(indx(ii)) !faintest galaxy
                rabs_n=rrmag(indx(4)) !faintest galaxy
                deltamag=rrmag(indx(ii))-rrmag(indx(1)) !faintest-brightest
                call indexx(ii,zc,indc)
                zc_max=zc(indc(ii))   !mayor redshift cosm
                zc_min=zc(indc(1))  !menor redshift cosm
                delta_zc=zc_max - zc_min 
                call indexx(ii,zs,inds)
                zs_max=zs(inds(ii))   !mayor redshift spect
                zs_min=zs(inds(1))  !menor redshift spect
                delta_zs=zs_max - zs_min 
                call indexx(ii,sstmass,indsat) !coeficiente satelites
                mst_1=sstmass(indsat(ii))  !1er stellar mass
                mst_2=sstmass(indsat(ii-1))  !2da stellar mass
                mst_3=sstmass(indsat(ii-2))  !3er stellar mass
                coef_mst=(mst_2+mst_3)/mst_1
                mst_1_bri=rrmag(indsat(ii))

                call indexx(ii,vvels,indv)
                vels_max=vvels(indv(ii))   !mayor vel miembro

                n3=0
                xxcm=0.
                yycm=0.
                zzcm=0.
                numx=0.
                denx=0.
                numy=0.
                deny=0.
                numz=0.
                denz=0.
                nd=0
                do k=1,nmie(j)
                 xx2=0.
                 !dist. entre glx
                 if(k.eq.nmie(j))goto 48  
                  do h=k+1,nmie(j)
                        nd=nd+1
                        xij=xx(k)-xx(h)
                        yij=yy(k)-yy(h)
                        zij=zz(k)-zz(h)
                        !modifico el "h"
                        if(xij.gt.lb2)xx(h)=xx(h)+lbox
                        if(yij.gt.lb2)yy(h)=yy(h)+lbox
                        if(zij.gt.lb2)zz(h)=zz(h)+lbox
                        if(xij.lt.-lb2)xx(h)=xx(h)-lbox
                        if(yij.lt.-lb2)yy(h)=yy(h)-lbox
                        if(zij.lt.-lb2)zz(h)=zz(h)-lbox
                        xij=xx(k)-xx(h)
                        yij=yy(k)-yy(h)
                        zij=zz(k)-zz(h)
                        dist_ij(nd)=sqrt(xij**2+yij**2+zij**2)
                        !xx2(k)=xx2(k) + xij !otra forma de clcular xcm

                  enddo
 48               continue 
                  if(rrmag(k).le.rabs1+3.)n3=n3+1  ! #miem en 3 mag


                  xxcm=xxcm+xx(k)
                  yycm=yycm+yy(k)
                  zzcm=zzcm+zz(k)

                  ! Ver donde va centro de masa (YOSE)
                  !pesados por masa ??
                  !xxcm=xxcm+xx(k)*masa(k)

C--------------- CALCULO XL DE LOS GRUPOS --------------------------
                  numx=numx+(xx(k)*lr(k))
                  denx=denx+lr(k)
                  numy=numy+(yy(k)*lr(k))
                  deny=deny+lr(k)
                  numz=numz+(zz(k)*lr(k))
                  denz=denz+lr(k)
c------------------------------------------------------------------

                  if(ncosmo.eq.3.or.ncosmo.eq.7.or.ncosmo.eq.8)then
                    write(90,26)ggID(k),xx(k),yy(k),zz(k),
     &                         sqrt(xx(k)**2+yy(k)**2+zz(k)**2),
     &                         vvx(k),vvy(k),vvz(k),
     &                         uumag(k),ggmag(k),rrmag(k),iimag(k),zzmag(k),
     &                         kkmag(k),jjmag(k),
     &                         sstmass(k),bbmass(k),iigru(k),iinum(k),k,
     &                         zc(k),zs(k),iigru_df(k),iinum_df(k)
                  else  
                    write(90,25)ggID(k),xx(k),yy(k),zz(k),
     &                         sqrt(xx(k)**2+yy(k)**2+zz(k)**2),
     &                         vvx(k),vvy(k),vvz(k),
     &                         uumag(k),ggmag(k),rrmag(k),iimag(k),zzmag(k),
     &                         sstmass(k),bbmass(k),iigru(k),iinum(k),k,
     &                         zc(k),zs(k),iigru_df(k),iinum_df(k)
                  endif
                end do

c-----------------------------------------------------------------------------------
                call mediana(dist_ij,nd,dist_mediana,sig_med) !dist. mediana entre glx
                dist_media=0.
                do i=1,nd
                 dist_media = dist_media + dist_ij(i)
                enddo
                dist_media=dist_media/real(nd) !dist. media entre glx

                !CM entre [-L/2, 3L/2]
                xcm=xxcm/real(nmie(j))
                ycm=yycm/real(nmie(j))
                zcm=zzcm/real(nmie(j))

                !otra forma de clcular xcm
                !xcm1=xx(1)-(xx2(1)/nmie(j))
                !write(177,*)j,xcm1

                xl=numx/denx !!!!! falta ver que dio
                yl=numy/deny !!!!! falta ver que dio
                zl=numz/denz !!!!! falta ver que dio


                call vecino(fac,j,ngrid,lbox,xcm,ycm,zcm,rvir_bi,rr_new,igru_new,lirstn,lgn,dmin,nvec,den_a,xvec,yvec,zvec)


                !densidad numerica local en cascaron (1<rvir_bi<3) sin miembros
                cascaron = 4.*pi*((3*rvir_bi)**3-(rvir_bi)**3)/3.
                den_aisl = den_a/cascaron

                den_cmp = nmie(j)/(4.*pi*((3*rvir_bi)))

                !write(106,*)j,den_cmp,den_aisl 

                
                !------------PERFILES------------------
                !perfil todas glx(miem, gal en gru , campo)
                call perfil2(0,fac,nbin,ngrid,rvir_bi,lbox,xcm,ycm,zcm,rr_new,stm_new,
     2          j,igru_new,lirstn,lgn,xbin,dens_new,dmax_x)
c
                !perfil glx miembro
                call perfil2(2,fac,nbin,ngrid,rvir_bi,lbox,xcm,ycm,zcm,rr_new,stm_new,
     2          j,igru_new,lirstn,lgn,xbin,dens_mie,dmax)
c
                !perfil glx campo
                call perfil2(1,fac,nbin,ngrid,rvir_bi,lbox,xcm,ycm,zcm,rr_new,stm_new,
     2          j,igru_new,lirstn,lgn,xbin,dens_cmp,dmax_x)

                !perfil glx otros grupos
                call perfil2(3,fac,nbin,ngrid,rvir_bi,lbox,xcm,ycm,zcm,rr_new,stm_new,
     2          j,igru_new,lirstn,lgn,xbin,dens_gru,dmax_x)



                write(99,77)j,xcm,ycm,zcm,vxcm,vycm,vzcm,
     &                      sig_3d,nmie(j),rvir_bi,mass_bi,
     &                      rabs1,rabs_n,deltamag,n3,
     &                      dist_media,dist_mediana,
     &                      nvec,rmag_new(nvec),xvec,
     &                      yvec,zvec,dmin,dmax,
     &                      zc_min,zc_max,delta_zc,           
     &                      zs_min,zs_max,delta_zs,
     &                      mst_1,mst_2,mst_3,mst_1_bri,
     &                      vels_max


                

                write(33,78)j,(xbin(k),k=1,nbin),(dens_cmp(k),k=1,nbin),(dens_new(k),k=1,nbin),
     & (dens_mie(k),k=1,nbin),(dens_gru(k),k=1,nbin)

        !dens_new TODO - (no)MIEMBROS (0,ROJO)        
        !dens_cmp CAMPO (1,AZUL)        
        !dens_mie MIEMBROS (2,NEGRO)     
        !dens_gru OTROS GRUPOS (3,VERDE)     

        if(mark(nvec).eq.1)nrep=nrep+1
        mark(nvec)=1
 12     end do

        write(*,*)'grupos sobredensos nmi>=3',ngfin

 23     format(i18,1x,7(f16.8,1x),5(f14.8,1x),2(e12.4,1x),4(i8,1x))
 24     format(i18,1x,7(f16.8,1x),7(f14.8,1x),2(e12.4,1x),4(i8,1x))

 25     format(i18,1x,7(f14.8,1x),5(f9.5,1x),2(e12.4,1x),2(i8,1x),i8,1x,2(f9.5,1x),2(i8,1x))
 26     format(i18,1x,7(f14.8,1x),7(f9.5,1x),2(e12.4,1x),2(i8,1x),i8,1x,2(f9.5,1x),2(i8,1x))

 77     format(i8,1x,3(f13.8,1x),3(f12.5,1x),f10.3,1x,i8
     &         ,1x,f7.4,1x,e12.4,1x,3(f9.5,1x),i8,1x
     & ,2(f10.5,1x),i9,1x,f17.8,5(1x,f10.5),6(1x,f9.5),3(1x,e12.4),
     &  1x,f9.5,1x,f13.5)

 78     format(i8,1x,75(f10.4,1x))
        
        close(99)
        close(90)
        close(33)
        end
        include '/big4/users/subs/fortran-doble/DLOCATE.FOR' 
        include '/big4/users/subs/fortran-doble/DINDEXX.FOR' 
        include '/big4/users/subs/porcentajes.f'
        include '/big4/users/subs/properties.f'
 
c-----------------------------------------------------------------------------
       subroutine radio_vir(ngal,xx,yy,zz,rvir,zmed)
       implicit none
       integer nmax,nn,ngal,i,l,k,indz
c       parameter(nmax=500000,nn=100000)
       real r(3,ngal),rcm(3,ngal),dis(ngal)
       real xx(ngal),yy(ngal),zz(ngal),rvir
       real pi,rrrr,rcc,x,y,z,zmed
       real*8 rmed
       real sumdis,disii,disjj,ctita,disd,stita,tita
       integer ntab
       real ztab(70000)
       real rtab(70000)
       common /tabla/ ztab,rtab,ntab  
c
       pi=4.*atan(1.)
c
       rmed=0.
       do i=1,ngal
        rmed=rmed+sqrt(xx(i)**2+yy(i)**2+zz(i)**2)
        r(1,i)=xx(i)
        r(2,i)=yy(i)
        r(3,i)=zz(i)
       end do
       rmed=rmed/float(ngal)
       call locate(rtab,ntab,rmed,indz)
       rrrr=(ztab(indz+1)-ztab(indz))/(rtab(indz+1)-rtab(indz))
       zmed=rrrr*(rmed-rtab(indz))+ztab(indz)
c
       sumdis=0.
       do i=1,ngal
        disii=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
        do l=1,ngal
            disjj=sqrt(r(1,l)**2+r(2,l)**2+r(3,l)**2)
            if(i.eq.l)goto 40
            ctita=0.
            disd=0.
            do k=1,3
             ctita=ctita+r(k,l)*r(k,i)
            end do
            ctita=ctita/(disii*disjj)
            if(ctita.ge.1.) then 
             stita=0.0
             tita=0.0
c             write(*,*)'111111111111111111'
            else
             stita=sqrt(1.0-ctita**2)
             tita=atan2(stita,ctita)
            end if
            if(tita.lt.0.) tita=tita+2.0*pi
            disd=2.0*rmed*tan(tita/2.0)
            sumdis=sumdis+1./(disd+.00001)
 40     end do
       end do
       sumdis=sumdis/2./(float(ngal)*(float(ngal)-1))
       rvir=1./sumdis*pi/2.
c
       return
       end
c-----------------------------------------------------------------------------
      subroutine rvir3d_old(ngal,lb,xx,yy,zz,rvir_3d) 
      implicit none
      integer ngal,i,j,nn
c      parameter (nn=100000,jj=500000)
      real xx(ngal),yy(ngal),zz(ngal),rvir_3d
      real difx,dify,difz,lb,r,suma

        suma=0.
        difx=0.
        dify=0.
        difz=0.
        r=0
        do i=1,ngal-1
           do j=i+1,ngal
                difx=xx(i)-xx(j)
                dify=yy(i)-yy(j)
                difz=zz(i)-zz(j)

                if(difx.lt.-lb/2.)difx=difx+lb
                if(dify.lt.-lb/2.)dify=dify+lb
                if(difz.lt.-lb/2.)difz=difz+lb

                if(difx.gt.lb/2.)difx=difx-lb
                if(dify.gt.lb/2.)dify=dify-lb
                if(difz.gt.lb/2.)difz=difz-lb

                r=sqrt(difx**2.+dify**2.+difz**2.)
                suma=suma+1./r
          end do
        end do
        rvir_3d=real(ngal*(ngal-1))*2./2./suma
        return
        end
c-----------------------------------------------------------------------------
      subroutine rvir3d(ngal,lb,xx,yy,zz,mass,rvir_3d) 
      !La marca decide que galaxias no van a formar parte del 
      !calculo del radio por estar demasiado cerca de otra galaxias. 
      !Elimina la galaxia menos masiva.      
      implicit none
      integer ngal,i,j,jj,no_masiva,nngal
      integer marca(ngal),marca_g
      real xx(ngal),yy(ngal),zz(ngal),rvir_3d
      real difx,dify,difz,lb,r,suma
      real mass(ngal),ratio

        suma=0.
        difx=0.
        dify=0.
        difz=0.
        marca=0
        marca_g=0
        jj=0
        do i=1,ngal-1
          if(marca(i).eq.1)goto 55
           do j=i+1,ngal
             if(marca(j).eq.1)goto 54
                difx=xx(i)-xx(j)
                dify=yy(i)-yy(j)
                difz=zz(i)-zz(j)

                if(difx.lt.-lb/2.)difx=difx+lb
                if(dify.lt.-lb/2.)dify=dify+lb
                if(difz.lt.-lb/2.)difz=difz+lb

                if(difx.gt.lb/2.)difx=difx-lb
                if(dify.gt.lb/2.)dify=dify-lb
                if(difz.gt.lb/2.)difz=difz-lb

                if(mass(i).gt.mass(j))then
                   ratio=mass(i)/mass(j)
                   no_masiva=j
                else
                   ratio=mass(j)/mass(i)
                   no_masiva=i
                endif   

                !print*, i,j,no_masiva
                r=sqrt(difx**2.+dify**2.+difz**2.)

                !if(difx.lt.1.e-2.and.ratio.gt.5.)then !no lo uso
                if(r.lt.1.e-2)then
                   jj=jj+1
                   marca(no_masiva)=1
                   goto 54
                endif   
                suma=suma+1./r
                !print*, r
 54       end do
 55     end do
        if((ngal-jj).lt.4) marca_g=1
        nngal=ngal-jj
        if(suma.ne.0.)rvir_3d=real(nngal*(nngal-1))*2./2./suma
        if(suma.eq.0.)rvir_3d=-99.        
        return
        end
c-----------------------------------------------------------------------------
        subroutine dist_rij(ngal,lb,xx,yy,zz,nn,rij)
        implicit none
        integer i,j,ngal,nn
        real lb
        real rij(ngal**2)
        real xx(ngal),yy(ngal),zz(ngal)
        real difx,dify,difz,r,suma
        real mass(ngal),ratio

        suma=0.
        difx=0.
        dify=0.
        difz=0.
        r=0
        nn=0
        do i=1,ngal-1
           do j=i+1,ngal
                nn=nn+1
                difx=xx(i)-xx(j)
                dify=yy(i)-yy(j)
                difz=zz(i)-zz(j)

                if(difx.lt.-lb/2.)difx=difx+lb
                if(dify.lt.-lb/2.)dify=dify+lb
                if(difz.lt.-lb/2.)difz=difz+lb

                if(difx.gt.lb/2.)difx=difx-lb
                if(dify.gt.lb/2.)dify=dify-lb
                if(difz.gt.lb/2.)difz=difz-lb

                r=sqrt(difx**2.+dify**2.+difz**2.)
                rij(nn)=r
        end do
      end do
      end          
c-----------------------------------------------------------------------------
        subroutine sig3d(ngal,vvx,vvy,vvz,vxcm,vycm,vzcm,sig_3d)
        implicit none
        integer ngal,i,nn
c        parameter (nn=100000)
        real vvx(ngal),vvy(ngal),vvz(ngal),vxcm,vycm,vzcm,modul(ngal)
        real sig_3d,suma

        suma=0
        do i=1,ngal
           modul(i)=(vvx(i)-vxcm)**2+(vvy(i)-vycm)**2+(vvz(i)-vzcm)**2
           suma=suma+modul(i)
        end do
        sig_3d=sqrt(suma/real(ngal-1))
        return
        end
c-----------------------------------------------------------------------
        function iargc()
c     get number of command line arguments using gfortran
      implicit none
c output
c     number of command line arguments
      integer iargc
c executable
       iargc=command_argument_count()
      return
      end
c-----------------------------------------------------------------------
        subroutine mediana(prop,n,med,sig)
        implicit none
        integer n,j
        real*8 prop(n)
        real   med,sig,x,q1,q3
        integer indx(n)
        integer k1,k2,k3,n1,n2

        call indexx(n,prop,indx)
       
        k2=int(real(n)/2.)
        if(mod(n,2).eq.0)then
           med=(prop(indx(k2))+prop(indx(k2+1)))/2.
        else
           med=prop(indx(k2+1))
        endif
        
        k1=int(real(n)/4.)
        if(mod(n,4).eq.0)then
           q1=(prop(indx(k1))+prop(indx(k1+1)))/2.
        else
           q1=prop(indx(k1+1))
        endif
  
        k3=int(real(n)*3./4.)
        n1=3*n
        n2=4

        if(mod(n1,n2).eq.0)then
           q3=(prop(indx(k3))+prop(indx(k3+1)))/2.
        else
           q3=prop(indx(k3+1))
        endif


        sig=1.58*(q3-q1)/sqrt(real(n))


        return
        end
c--------------------------------------------------------------------
      subroutine grid(r,ngrid,nobj,lirst,lg)
      integer ngridmax,nmax
      parameter (ngridmax=128,nmax=15000000)
      integer lirst(ngridmax,ngridmax,ngridmax),lg(nmax)
      real r(3,nmax)
c
      write(*,*) 'grid: start'
c
      lirst=0
      lg=0
c
      do i=1,nobj
        indx=int(r(1,i)+1.)
        indy=int(r(2,i)+1.)
        indz=int(r(3,i)+1.)
        lirst(indx,indy,indz)=i
      end do
c
      do i=1,nobj
        indx=int(r(1,i)+1.)
        indy=int(r(2,i)+1.)
        indz=int(r(3,i)+1.)
        lg(lirst(indx,indy,indz))=i
        lirst(indx,indy,indz)=i
      end do
c
      write(*,*) 'grid: end'
      return
      end
c---------------------------------------
        subroutine vecino(fac,ig,ngrid,lbox,xcm,ycm,zcm,rvir,rr,igru,lirst,lg,dmin,nvec,den,xvec,yvec,zvec)
        implicit none
        integer nmax,ngrid,nngrid
        parameter (nmax=15000000,nngrid=128)
        integer lirst(nngrid,nngrid,nngrid),lg(nmax)
        real fac,xcm,ycm,zcm,lbox,lbox2,pi
        real xx,yy,zz,xc,yc,zc,xii,yii,zii
        integer ipx,ipy,ipz,ix,iy,iz,ii,ixt,iyt,izt
        real rr(3,nmax),rvir
        integer ig,igru(nmax),nvec
        real dd,dmin,den
        real xvec,yvec,zvec


        lbox2=lbox/2.
        dmin=1.e10
        xx=fac*xcm
        yy=fac*ycm
        zz=fac*zcm
        ipx=int(xx+1.)
        ipy=int(yy+1.)
        ipz=int(zz+1.)
c---------------------------------
         do ix=ipx-1,ipx+1
           if(ix.eq.0)then
                ixt=ngrid
           else if(ix.eq.ngrid+1)then
                ixt=1
           else
                ixt=ix
           endif
          do iy=ipy-1,ipy+1
            if(iy.eq.0)then
                iyt=ngrid
            else if(iy.eq.ngrid+1)then
                iyt=1
            else
                iyt=iy
            endif
            do iz=ipz-1,ipz+1
              if(iz.eq.0)then
                izt=ngrid
              else if(iz.eq.ngrid+1)then
                izt=1
              else
                izt=iz
              endif

              den=0. 
              ii=lirst(ixt,iyt,izt)
              if(ii.eq.0) goto 31
 21           ii=lg(ii)
              if(igru(ii).eq.ig)goto 15 !elimina los miembros 
              xii=rr(1,ii)/fac 
              yii=rr(2,ii)/fac 
              zii=rr(3,ii)/fac 
              if((xcm-xii).gt.lbox2)xii=xii+lbox
              if((ycm-yii).gt.lbox2)yii=yii+lbox
              if((zcm-zii).gt.lbox2)zii=zii+lbox
              if((xcm-xii).lt.-lbox2)xii=xii-lbox
              if((ycm-yii).lt.-lbox2)yii=yii-lbox
              if((zcm-zii).lt.-lbox2)zii=zii-lbox
              dd=sqrt((xcm-xii)**2.+(ycm-yii)**2.+(zcm-zii)**2.)
              if(dd.lt.dmin)then
                dmin=dd
                xvec=xii
                yvec=yii
                zvec=zii
                nvec=ii
              endif
c-------      !densidad local en cascaron (1<rvir_3d<3)
              if(dd.gt.rvir/fac.and.dd.lt.3*rvir/fac) den=den+1 

 15           if(ii.eq.lirst(ixt,iyt,izt))goto 31
              goto 21
c-------
c 15           if(ii.eq.lirst(ixt,iyt,izt)) goto 31
c              goto 21
 31         end do
          end do
         end do
        

        return
        end

c---------------------------------------
        subroutine perfil(fac,nbin,ngrid,rvir,lbox,xcm,ycm,zcm,rr,
     1  lirst,lg,xbin,dens)
        implicit none
        integer nmax,ngrid,nngrid,nbin,i
        parameter (nmax=15000000,nngrid=128)
        integer lirst(nngrid,nngrid,nngrid),lg(nmax)
        real fac,xcm,ycm,zcm,lbox,lbox2,rvir,pi
        real xx,yy,zz,xc,yc,zc,xii,yii,zii
        integer ipx,ipy,ipz,ix,iy,iz,ii,ixt,iyt,izt
        real rr(3,nmax)
        integer ig,igru(nmax),nvec
        real dd,dmin,dbin
        integer ibin,nb(nbin)
        real xmin,xmax,vol(nbin)
        real xbin(nbin),dens(nbin)

        pi=4.*atan(1.)

        dbin=5./real(nbin) 

        do i=1,nbin
           xbin(i)=(real(i)-0.5)*dbin
           xmin=(real(i)-1.)*dbin
           xmax=real(i)*dbin
           vol(i)=4.*pi/3.*(xmax**3-xmin**3)
        end do

        
        nb=0
        lbox2=lbox/2.
        dmin=1.e10
        xx=fac*xcm
        yy=fac*ycm
        zz=fac*zcm
        ipx=int(xx+1.)
        ipy=int(yy+1.)
        ipz=int(zz+1.)
c---------------------------------
         do ix=ipx-1,ipx+1
           if(ix.eq.0)then
                ixt=ngrid
           else if(ix.eq.ngrid+1)then
                ixt=1
           else
                ixt=ix
           endif
          do iy=ipy-1,ipy+1
            if(iy.eq.0)then
                iyt=ngrid
            else if(iy.eq.ngrid+1)then
                iyt=1
            else
                iyt=iy
            endif
            do iz=ipz-1,ipz+1
              if(iz.eq.0)then
                izt=ngrid
              else if(iz.eq.ngrid+1)then
                izt=1
              else
                izt=iz
              endif

              ii=lirst(ixt,iyt,izt)
              if(ii.eq.0) goto 31
 21           ii=lg(ii)
c              if(igru(ii).eq.ig)goto 15 !elimina los miembros 
              xii=rr(1,ii)/fac 
              yii=rr(2,ii)/fac 
              zii=rr(3,ii)/fac 
              if((xcm-xii).gt.lbox2)xii=xii+lbox
              if((ycm-yii).gt.lbox2)yii=yii+lbox
              if((zcm-zii).gt.lbox2)zii=zii+lbox
              if((xcm-xii).lt.-lbox2)xii=xii-lbox
              if((ycm-yii).lt.-lbox2)yii=yii-lbox
              if((zcm-zii).lt.-lbox2)zii=zii-lbox
              dd=sqrt((xcm-xii)**2.+(ycm-yii)**2.+(zcm-zii)**2.)/rvir
              if(dd.gt.5.)goto 15
              ibin=int(dd/dbin)+1
              if(ibin.eq.nbin+1)ibin=nbin
              nb(ibin)=nb(ibin)+1
              
 15           if(ii.eq.lirst(ixt,iyt,izt)) goto 31
              goto 21
 31         end do
          end do
         end do
        do i=1,nbin
           dens(i)=real(nb(i))/vol(i)
        end do


        return
        end

c        ---------------------------------------
        subroutine perfil2(ind,fac,nbin,ngrid,rvir,lbox,xcm,ycm,zcm,rr,stm,
     1  ig,igru,lirst,lg,xbin,dens,dmax)
        implicit none
        integer nmax,ngrid,nngrid,nbin,i
        parameter (nmax=15000000,nngrid=128)
        integer lirst(nngrid,nngrid,nngrid),lg(nmax)
        real fac,xcm,ycm,zcm,lbox,lbox2,rvir,pi
        real xx,yy,zz,xc,yc,zc,xii,yii,zii
        integer ipx,ipy,ipz,ix,iy,iz,ii,ixt,iyt,izt
        real rr(3,nmax),stm(nmax)
        integer ig,igru(nmax),nvec,ind
        real dd,dmin,dbin
        integer ibin,nb(nbin)
        real rnb(nbin)
        real xmin,xmax,vol(nbin)
        real xbin(nbin),dens(nbin)
        real dmax

        pi=4.*atan(1.)

        dbin=5./real(nbin) 

        do i=1,nbin
           xbin(i)=(real(i)-0.5)*dbin
           xmin=(real(i)-1.)*dbin
           xmax=real(i)*dbin
           vol(i)=4.*pi/3.*(xmax**3-xmin**3)
        end do

        dmax=0.
        nb=0
        rnb=0.
        lbox2=lbox/2.
        dmin=1.e10
        xx=fac*xcm
        yy=fac*ycm
        zz=fac*zcm
        ipx=int(xx+1.)
        ipy=int(yy+1.)
        ipz=int(zz+1.)
c---------------------------------
         do ix=ipx-1,ipx+1
           if(ix.eq.0)then
                ixt=ngrid
           else if(ix.eq.ngrid+1)then
                ixt=1
           else
                ixt=ix
           endif
          do iy=ipy-1,ipy+1
            if(iy.eq.0)then
                iyt=ngrid
            else if(iy.eq.ngrid+1)then
                iyt=1
            else
                iyt=iy
            endif
            do iz=ipz-1,ipz+1
              if(iz.eq.0)then
                izt=ngrid
              else if(iz.eq.ngrid+1)then
                izt=1
              else
                izt=iz
              endif

              ii=lirst(ixt,iyt,izt)
              if(ii.eq.0) goto 31
 21           ii=lg(ii)
              if(ind.eq.0.and.igru(ii).eq.ig)goto 15 !todo menos miembros
              if(ind.eq.1.and.igru(ii).ne.0)goto 15  !solo campo 
              if(ind.eq.2.and.igru(ii).ne.ig)goto 15 !solo miembros
              if(ind.eq.3.and.igru(ii).eq.ig)goto 15 !solo gal in gru
              if(ind.eq.3.and.igru(ii).eq.0)goto 15  !solo gal in gru
              xii=rr(1,ii)/fac 
              yii=rr(2,ii)/fac 
              zii=rr(3,ii)/fac 
              if((xcm-xii).gt.lbox2)xii=xii+lbox
              if((ycm-yii).gt.lbox2)yii=yii+lbox
              if((zcm-zii).gt.lbox2)zii=zii+lbox
              if((xcm-xii).lt.-lbox2)xii=xii-lbox
              if((ycm-yii).lt.-lbox2)yii=yii-lbox
              if((zcm-zii).lt.-lbox2)zii=zii-lbox
              dd=sqrt((xcm-xii)**2.+(ycm-yii)**2.+(zcm-zii)**2.)/rvir
              if(ind.eq.2)then
                 if(dd.gt.dmax)dmax=dd
              endif
              if(dd.gt.5.)goto 15
              ibin=int(dd/dbin)+1
              if(ibin.eq.nbin+1)ibin=nbin
              nb(ibin)=nb(ibin)+1
              !rnb(ibin)=rnb(ibin)+stm(ii)
              
 15           if(ii.eq.lirst(ixt,iyt,izt)) goto 31
              goto 21
 31         end do
          end do
         end do
        do i=1,nbin
           dens(i)=real(nb(i))/vol(i)
           !dens(i)=rnb(i)/vol(i)
        end do


        return
        end

c---------------------------------------
        subroutine den_stm(ind,fac,nbin,ngrid,rvir,lbox,xcm,ycm,zcm,rr,stm,
     1  ig,igru,lirst,lg,xbin,dens,dmax)
        implicit none
        integer nmax,ngrid,nngrid,nbin,i
        parameter (nmax=15000000,nngrid=128)
        integer lirst(nngrid,nngrid,nngrid),lg(nmax)
        real fac,xcm,ycm,zcm,lbox,lbox2,rvir,pi
        real xx,yy,zz,xc,yc,zc,xii,yii,zii
        integer ipx,ipy,ipz,ix,iy,iz,ii,ixt,iyt,izt
        real rr(3,nmax),stm(nmax)
        integer ig,igru(nmax),nvec,ind
        real dd,dmin,dbin
        integer ibin,nb(nbin)
        real rnb(nbin)
        real xmin,xmax,vol(nbin)
        real xbin(nbin),dens(nbin)
        real dmax

        pi=4.*atan(1.)

        dbin=5./real(nbin) 

        do i=1,nbin
           xbin(i)=(real(i)-0.5)*dbin
           xmin=(real(i)-1.)*dbin
           xmax=real(i)*dbin
           vol(i)=4.*pi/3.*(xmax**3-xmin**3)
        end do

        dmax=0.
        nb=0
        rnb=0.
        lbox2=lbox/2.
        dmin=1.e10
        xx=fac*xcm
        yy=fac*ycm
        zz=fac*zcm
        ipx=int(xx+1.)
        ipy=int(yy+1.)
        ipz=int(zz+1.)
c---------------------------------
         do ix=ipx-1,ipx+1
           if(ix.eq.0)then
                ixt=ngrid
           else if(ix.eq.ngrid+1)then
                ixt=1
           else
                ixt=ix
           endif
          do iy=ipy-1,ipy+1
            if(iy.eq.0)then
                iyt=ngrid
            else if(iy.eq.ngrid+1)then
                iyt=1
            else
                iyt=iy
            endif
            do iz=ipz-1,ipz+1
              if(iz.eq.0)then
                izt=ngrid
              else if(iz.eq.ngrid+1)then
                izt=1
              else
                izt=iz
              endif

              ii=lirst(ixt,iyt,izt)
              if(ii.eq.0) goto 31
 21           ii=lg(ii)
              if(ind.eq.0.and.igru(ii).eq.ig)goto 15 !todo menos miembros
              if(ind.eq.1.and.igru(ii).ne.0)goto 15  !solo campo 
              if(ind.eq.2.and.igru(ii).ne.ig)goto 15 !solo miembros
              if(ind.eq.3.and.igru(ii).eq.ig)goto 15 !solo gal in gru
              if(ind.eq.3.and.igru(ii).eq.0)goto 15  !solo gal in gru
              xii=rr(1,ii)/fac 
              yii=rr(2,ii)/fac 
              zii=rr(3,ii)/fac 
              if((xcm-xii).gt.lbox2)xii=xii+lbox
              if((ycm-yii).gt.lbox2)yii=yii+lbox
              if((zcm-zii).gt.lbox2)zii=zii+lbox
              if((xcm-xii).lt.-lbox2)xii=xii-lbox
              if((ycm-yii).lt.-lbox2)yii=yii-lbox
              if((zcm-zii).lt.-lbox2)zii=zii-lbox
              dd=sqrt((xcm-xii)**2.+(ycm-yii)**2.+(zcm-zii)**2.)/rvir
              if(ind.eq.2)then
                 if(dd.gt.dmax)dmax=dd
              endif
              if(dd.gt.5.)goto 15
              ibin=int(dd/dbin)+1
              if(ibin.eq.nbin+1)ibin=nbin
              nb(ibin)=nb(ibin)+1
              !rnb(ibin)=rnb(ibin)+stm(ii)
              
 15           if(ii.eq.lirst(ixt,iyt,izt)) goto 31
              goto 21
 31         end do
          end do
         end do
        do i=1,nbin
           dens(i)=real(nb(i))/vol(i)
           !dens(i)=rnb(i)/vol(i)
        end do


        return
        end


