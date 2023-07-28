c------gfortran -o mock_difusos mock_difusos.f
c input1: cosmo
       !---------------------------------
       module files
       real lbox,h0
       character*79 filein1,filein2,filein3,filecosmo
       character*79 fileout1,fileout2,fileout3,fileout4
       end module files
       !---------------------------------
       program mock_ais
       use files        
       implicit none 
       integer i,j,h,nlines,nlim
       integer nmax_1,ngru,ngal,nmax_gal,igru_lim,igru_lim_2
       parameter (nmax_1=90000000) !rango igru(boxes) [10000059,81023301]
       parameter (nmax_gal=1000000000) !numero de galaxias (boxes)
       parameter (igru_lim= 10000000) !10 millones (base de igru)
       parameter (igru_lim_2= 100000) !100000      (base de igru MII)
       integer igrumax,ngalbox,igrumax8
       real den_num_med
       !rango gID(box) [69954,510000127000295]


       !grupos
       integer j1  !nombre del grupo
       real xxcm,yycm,zzcm
       real vxcm,vycm,vzcm
       real xcm(nmax_1),ycm(nmax_1),zcm(nmax_1)
       !integer nmie,n3
       !real sig_3d,rvir,mvir_bi
       !real rabs1,rabs_n,deltamag
       integer,allocatable,dimension(:) :: mark
       integer marca(nmax_1)
       integer nm(nmax_1),nmi(nmax_1)
       integer*8 idgru
       character*86 resto(nmax_1),resto1
       real alfacm(nmax_1),delcm(nmax_1),zspec_cm(nmax_1)
       real rproycm,rocm,zcosmo,vpcm

       !galaxias
       integer*8 gIDt   !ID de la galaxia en el box
       integer*8 galID  !nombre de las galaxias cambiado p/ 8 boxes
       real x,y,z,dist
       real difx,dify,difz
       real vx,vy,vz
       real umag,gmag,rmag,imag,zmag
       real kmag,jmag
       real stmass,bmass
       integer inum,k
       integer igru

       integer indr,indz
       real zc,zs,zzz,vp,vr
       real s,rrr
       real alfa,del,ro,rproy
       real c,pi
       real kcor,kcorr,color_gr
       real mlim,r_ap,g_ap

       integer ng1,ng2,ng3,ng4,nmmax,nmimax,n10,nz0
       integer nn,nmock,ngal_ais

       real kcorru,kcorrg,kcorrr,kcorri,kcorrz
       real gobs,robs

       !lazo box
       integer bx,by,bz,jj
       !real lbox,h0
       integer nbox,iigru
       parameter (nbox=1)
       real xx,yy,zz


c      Tabla redshift
       integer ntab
       real ztab(70000)
       real rtab(70000)
       !common /tabla/ ztab,rtab,ntab

c-----------------------------------------
       external  iargc
       integer iargc,snap,ncosmo
       character*80 string

      !----------------------------------------
      !----------------------------------------
       c = 299792.458 !Km/s
       pi = 4.*atan(1.)
       mlim = 17.77

       if (iargc().ne.1) then
       print *,'usage:  ncosmo '
       print *,'Wmap_1/guo11 1 '
       print *,'Wmap_7/guo13 2 '
       print *,'Plank/hen15 3 '
       print *,'Wmap_1/guo11-MII 4 '
       stop
       endif

       call getarg(1,string)
       read(string,*)ncosmo 
  
       call reader_files(ncosmo) 
       open(30,file=filein1,status='old')
       open(31,file=filein2,status='old')
       open(33,file=filein3,status='old')
       open(40,file=fileout1,status='unknown')
       open(89,file=filecosmo,status='old')
       open(90,file=fileout3,status='unknown')
       open(91,file=fileout4,status='unknown')


      !----------------------------------------
       read(33,*)igrumax,ngalbox,den_num_med
       if((igru_lim+igrumax).gt.nmax_1)then
           print*,'mal dimensiones!!!!!!!!'
           stop
       endif    
      !----------------------------------------
      !----------------------------------------
       call count_lines(89,nlines)
       ntab = nlines
       print*,'tabla de dsit',ntab
       do i=1,ntab
            read(89,*)ztab(i),rtab(i)
       end do
       close(89)
      !----------------------------------------
      !----------------------------------------
       allocate(mark(igrumax))    !rango iigru(box)  [59,1023301]

       call count_lines(30,nlines)
       ngru = nlines

       mark = 0  !1box
       marca = 0 !8boxes
       igrumax8=0
       do i=1,ngru
            read(30,'(i8,1x,3(f10.5,1x),3(f12.5,1x),a86)')
     &                                  j1,xxcm,yycm,zzcm,
     &                                  vxcm,vycm,vzcm,resto1

            !character resto1 estan las propiedades 3d
            ! sig_3d,nmie,rvir,mvir_bi,rabs1,rabs_n,deltamag,n3
            

            mark(j1)= 1

            !pego boxes
            do bx=1,nbox
             do by=1,nbox
              do bz=1,nbox

                jj=bx+(by-1)*nbox+(bz-1)*nbox*nbox !nombres(#) boxes
                if(ncosmo.eq.4)then
                       idgru = j1 + igru_lim_2*jj
                else
                       idgru = j1 + igru_lim*jj
                endif

                if(idgru.GT.igrumax8)igrumax8=idgru


                xcm(idgru) = xxcm + real(bx-1)*lbox
                ycm(idgru) = yycm + real(by-1)*lbox
                zcm(idgru) = zzcm + real(bz-1)*lbox

                resto(idgru)=resto1

                !solo me quedo con los CG asilados en un octante
                !esferico
                rocm = sqrt(xcm(idgru)**2+ycm(idgru)**2+zcm(idgru)**2)  
                if(rocm.le.nbox*lbox)marca(idgru)=1

                !con  x, y, z calculo alfa y delta (pi/2-tita)
                rproycm = sqrt(xcm(idgru)**2+ycm(idgru)**2)
                delcm(idgru) = atan2(zcm(idgru),rproycm)
                alfacm(idgru) = atan2(ycm(idgru),xcm(idgru))
                !estoy en un octante
                if(alfacm(idgru).lt.0.)alfacm(idgru)=alfacm(idgru)+2.*pi

                !redshift cosmologico    
                call locate(rtab,ntab,rocm,indr)
                zzz=(ztab(indr+1)-ztab(indr))/(rtab(indr+1)-rtab(indr))
                zcosmo=zzz*(rocm-rtab(indr))+ztab(indr)  

                !velocidades pec=r.v/|r| y redshift cosmologico
                vpcm=(vxcm*xcm(idgru)+vycm*ycm(idgru)+vzcm*zcm(idgru))/rocm
                zspec_cm(idgru)=(1.+zcosmo)*(1.+vpcm/c)-1. 

              enddo   
             enddo   
            enddo   


       enddo
       close(30)
      !----------------------------------------
      !----------------------------------------
       call count_lines(31,nlines)
       ngal = nlines

       ngal_ais = 0
       nmock = 0
       nlim = 0
       nmi = 0
       nm = 0
       n10 = 0
       nz0 = 0
       do i=1,ngal
          !------------
          if(ncosmo.eq.1.or.ncosmo.eq.2.or.ncosmo.eq.4)then
            read(31,*)gIDt,xx,yy,zz,dist,vx,vy,vz,
     &                    umag,gmag,rmag,imag,zmag,
     &                    stmass,bmass,iigru,inum,k,
     &                    zc,zs
          elseif(ncosmo.eq.3)then
            read(31,*)gIDt,xx,yy,zz,dist,vx,vy,vz,
     &                    umag,gmag,rmag,imag,zmag,
     &                    kmag,jmag,
     &                    stmass,bmass,iigru,inum,k,
     &                    zc,zs
          endif
          !------------
   
            if(mark(iigru).eq.0) goto 57  !solo glx en difusos
            ngal_ais = ngal_ais + 1
            
            !pego boxes
            do bx=1,nbox
             do by=1,nbox
              do bz=1,nbox

                !buscar un buen nombre para los grupos
                !el igru_max(grupos aislados) = 1023301
                !le sumo 2 ordenes mas asi no se repiten los id
                jj=bx+(by-1)*nbox+(bz-1)*nbox*nbox !nombres(#) boxes
                if(ncosmo.eq.4)then
                       igru = iigru + igru_lim_2*jj
                else
                       igru = iigru + igru_lim*jj
                endif


                !¿necesito cambiar el gID?
                !write(85,*)gIDt
                galID=nmax_gal*jj + i
                

                if(marca(igru).eq.0)goto 58  !solo glx en aisl en un oct esferico
                nmock = nmock + 1
               

                x = xx + real(bx-1)*lbox
                y = yy + real(by-1)*lbox
                z = zz + real(bz-1)*lbox

                ro = sqrt(x**2+y**2+z**2)
                !write(85,*)ro                

                !con  x, y, z calculo alfa y delta (pi/2-tita)
                rproy = sqrt(x**2+y**2)
                del = atan2(z,rproy)
                alfa = atan2(y,x)
                !estoy en un octante
                if(alfa.lt.0.)alfa = alfa + 2.*pi

                !redshift cosmologico    
                call locate(rtab,ntab,ro,indr)
                zzz=(ztab(indr+1)-ztab(indr))/(rtab(indr+1)-rtab(indr))
                zc=zzz*(ro-rtab(indr))+ztab(indr)  

                !velocidades pec=r.v/|r| y redshift cosmologico
                vp=(vx*x+vy*y+vz*z)/ro
                zs=(1.+zc)*(1.+vp/c)-1. 

                if(zs.lt.0.)nz0=nz0+1
                if(zs.lt.0.)goto 58

                !velocidad radial
                vr=zs*c

                !distancia comovil distor    
                call locate(ztab,ntab,zs,indz)
                rrr=(rtab(indz+1)-rtab(indz))/(ztab(indz+1)-ztab(indz))
                s=rrr*(zs-ztab(indz))+rtab(indz)
    
                !empiezo el calculo de las Kcorr
                r_ap = rmag + 25. + 5.*log10(s*(1. + zs)) ! mag ap rest frame
                g_ap = gmag + 25. + 5.*log10(s*(1. + zs)) ! mag ap rest frame

                !chilinga (calculo K-corrction de forma iterativa)
                !output: magnitud aparente obs. frame (robs,gobs)
                call kcorr_it(2,3,g_ap,r_ap,zs,kcorrg,kcorrr,gobs,robs)

!                if(robs.lt.0.)goto 58
                if(robs.lt.10.)n10=n10+1

                nmi(igru) = nmi(igru) + 1 !nmi:numero de miembros

                if(robs.gt.mlim)goto 58  !solo glx en aisl en un oct esferico c/ corte mlim
                nlim = nlim + 1

                nm(igru) = nm(igru) + 1 !nm:numero de miembros con corte en mag


                 !------------
                 if(ncosmo.eq.1.or.ncosmo.eq.2.or.ncosmo.eq.4)then
                    write(40,77)galID,igru,inum,k,robs,gobs,
     &                          alfa*180./pi,del*180./pi,zs,
     &                          x,y,z,dist,vx,vy,vz,
     &                          umag,gmag,rmag,imag,zmag,
     &                          stmass,bmass

                 elseif(ncosmo.eq.3)then
                    write(40,78)galID,igru,inum,k,robs,gobs,
     &                          alfa*180./pi,del*180./pi,zs,
     &                          x,y,z,dist,vx,vy,vz,
     &                          umag,gmag,rmag,imag,zmag,
     &                          kmag,jmag,
     &                          stmass,bmass

                 endif  
                 !------------


 58           enddo
             enddo
            enddo !lazo boxes


 57    enddo  !lazo glx
       close(31)
       close(40)
       !close(41)

       nn=0
       ng1 = 0
       ng2 = 0
       ng3 = 0
       ng4 = 0
       nmmax = 0
       nmimax = 0

       do i =1,nmax_1
          if(marca(i).eq.1)then
              nn=nn+1
              write(90,79)i,xcm(i),ycm(i),zcm(i),
     &                    alfacm(i),delcm(i),zspec_cm(i),resto(i),nm(i)
          endif 

          if(nm(i).eq.1)ng1 = ng1 +1
          if(nm(i).eq.2)ng2 = ng2 +1
          if(nm(i).eq.3)ng3 = ng3 +1
          if(nm(i).gt.3)then
              ng4 = ng4 +1
              !grupos aislados c/corte mlim >4 miemb
              write(91,79)i,xcm(i),ycm(i),zcm(i),
     &                    alfacm(i),delcm(i),zspec_cm(i),resto(i),nm(i)
          endif
          if(nm(i).gt.nmmax)nmmax = nm(i)
          if(nmi(i).gt.nmimax)nmimax = nmi(i)
       enddo
       close(90)
       close(91)
       print*, '-------------------------'
       print*,'grupos aislados box',ngru
       print*, 'grupos_aisl octante de r 2*lbox',nn
       print*, 'grupos_aisl octante con corte en mag y nmi>4',ng4
       print*, 'grupos_aisl octante c/ corte en mag 1, 2 y 3 miembros',ng1,ng2,ng3
       print*, '-------------------------'
       print*, 'glx en grupos_aisl box',ngal_ais
       print*, 'glx en grupos_aisl octante de r 2*lbox',nmock
       print*, 'glx en grupos_aisl octante de r 2*lbox (m_ap < m_lim)',nlim
       print*, '-------------------------'
       print*,'glx en grupos_aisl con zs<0',nz0
       print*,'glx en grupos_aisl con 0<robs<10',n10
       print*, '-------------------------'
       print*, 'max numero de miembros',nmimax
       print*, 'max numero de miembros con corte en mag ',nmmax
       print*, '-------------------------'
       !print*,'galaxias en grupos sobredensos box',ngal
       print*, 'nombre maximo de grupos sobredensos 8 boxes',igrumax8

      !----------------------------------------
      !----------------------------------------
   
 77   format(i18,1x,3(i8,1x),2(f14.6,1x),1x,2(f14.8,1x),f14.7,1x,
     &    7(f14.6,1x),5(f14.6,1x),2(e12.4,1x))

 78   format(i18,1x,3(i8,1x),2(f14.6,1x),1x,2(f14.8,1x),f14.7,1x,
     &    7(f14.6,1x),7(f14.6,1x),2(e12.4,1x))


 79    format(i8,1x,3(f10.5,1x),3(f14.7,1x),a86,1x,i3)


       end program


      include 'kcorr_iterative.f' !(contine tmb al "locate")
c-----------------------------------------------------------------------
c-------------------SUBRUTINAS Y FUNCIONES------------------------------
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
       subroutine count_lines(nunit,nlines)
       implicit none
       integer, intent(in)            :: nunit
       integer, intent(out)           :: nlines
       integer                        :: check
        nlines = 0
        do
        read(nunit,*,iostat=check)
        if(check /= 0)exit
         nlines = nlines + 1
        end do
        rewind(nunit)
       end subroutine count_lines

c-----------------------------------------
       module constantes !ver
       implicit none
       real c

       !c=299792.458 !Km/s
       !G=6.6726*10.**(-11.)  !m**3/Kg s**2
       !munit=(3.0857*10.**(22.))/(1.99*10.**(30.))
       !lsun=3.846*(10.**33.) !erg/s
       !msun=4.65 !r-band http://mips.as.arizona.edu/~cnaw/sun.html 
       !pi=4.*atan(1.)
       end module constantes

c-----------------------------------------
       subroutine reader_files(ncosmo)
       use files
       implicit none
       integer ncosmo

       if(ncosmo.eq.1)then
         lbox=500.
         h0=0.73
         filecosmo='../Mill_I/red2dis_lcdm.dat-mill_wmap1'
           write(filein1,'("guo_11/tablagru.dat")')
           write(filein2,'("guo_11/galx_gru.dat")')
           write(filein3,'("guo_11/file_dim.dat")')
           write(fileout1,'("guo_11/gal_gru_mock.dat")')
           write(fileout3,'("guo_11/gru_mock.dat")')
           write(fileout4,'("guo_11/gru_mock_mlim.dat")')
       endif
       !---------------- 
       if(ncosmo.eq.2)then
         lbox=500.
         h0=0.704
         filecosmo='../Mill_I/red2dis_lcdm.dat-mill_wmap7'
           write(filein1,'("guo_13/tablagru.dat")')
           write(filein2,'("guo_13/galx_gru.dat")')
           write(filein3,'("guo_13/file_dim.dat")')
           write(fileout1,'("guo_13/gal_gru_mock.dat")')
           write(fileout3,'("guo_13/gru_mock.dat")')
           write(fileout4,'("guo_13/gru_mock_mlim.dat")')
       endif
       !----------------
       if(ncosmo.eq.3)then
         lbox=480.279
         h0=0.673
         filecosmo='../Mill_I/red2dis_lcdm.dat-mill_planck'
           write(filein1,'("hen_15/tablagru.dat")')
           write(filein2,'("hen_15/galx_gru.dat")')
           write(filein3,'("hen_15/file_dim.dat")')
           write(fileout1,'("hen_15/gal_gru_mock.dat")')
           write(fileout3,'("hen_15/gru_mock.dat")')
           write(fileout4,'("hen_15/gru_mock_mlim.dat")')
       endif
       !----------------
       if(ncosmo.eq.4)then
         lbox=100.
         h0=0.73
         filecosmo='../Mill_II/red2dis_lcdm.dat-mill_wmap1'
           write(filein1,'("guo_II/tablagru.dat")')
           write(filein2,'("guo_II/galx_gru.dat")')
           write(filein3,'("guo_II/file_dim.dat")')
           write(fileout1,'("guo_II/gal_gru_mock.dat")')
           write(fileout3,'("guo_II/gru_mock.dat")')
           write(fileout4,'("guo_II/gru_mock_mlim.dat")')
       endif
       !---------------- 

       end subroutine


