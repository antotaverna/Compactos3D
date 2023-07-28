!gfortran -o readgroups2 readgroups2.f
!input1: cosmologia 1, 2, 3 o 4
!./readgroups2 input1!

!El archivo out.dat de Hen_15 lee 2 lineas mas (kmag,jmag)
       program readdatos

        implicit none
        integer i,nrep,col1,col2,nmax,nglx
        parameter(nmax=15400000)
        integer rep(nmax),ifila,igrut,inumt,nfof
        integer igru(nmax),inum(nmax),ty

        integer*8 gID
        real x,y,z,vx,vy,vz,umag,gmag,rmag,imag,zmag,kmag,jmag
        real stmass,bmass,sfr
        real age,ubmag,gbmag,rbmag,ibmag,zbmag,dist

        integer snap,ncosmo
        integer t1,t2

        external  iargc
        integer iargc
        character*80 string
        
        character*80 fileglx,filerep,fileout,filefof
c   
        integer ipcall,iplast,iphead,jc
        common /kpercent/ ipcall,iplast,iphead,jc

        if (iargc().ne.1) then
        print *,'usage:  ncosmo'
        print *,'Wmap_1/guo11 1 '
        print *,'Wmap_7/guo13 2 '
        print *,'Plank/hen15 3 '        
        print *,'Wmap_1/guo11-MII 4 '
        print *,'Wmap_1/guo11-MII_16 5 '
        stop
        endif

        call getarg(1,string)
        read(string,*)ncosmo 

c------------ FILENAMES -------------------------------        
        if(ncosmo.eq.1)then
        write(fileglx,'("guo_11/tablaglx.dat")')
        write(filerep,'("../Mill_I/Guo11_wmap1/z0/repetidas09.dat-0")')
        write(fileout,'("../Mill_I/Guo11_wmap1/z0/out09.dat-0")')
        write(filefof,'("../Mill_I/Guo11_wmap1/z0/fof09.dat")')
        endif
        !----------------
        if(ncosmo.eq.2)then
        write(fileglx,'("guo_13/tablaglx.dat")')
        write(filerep,'("../Mill_I/Guo13_wmap7/z0/repetidas09.dat-0")')
        write(fileout,'("../Mill_I/Guo13_wmap7/z0/out09.dat-0")')
        write(filefof,'("../Mill_I/Guo13_wmap7/z0/fof09.dat")')
        endif
        !----------------
        if(ncosmo.eq.3)then
        write(fileglx,'("hen_15/tablaglx.dat")')
        write(filerep,'("../Mill_I/Hen15_plank/z0/repetidas09.dat-0")')
        write(fileout,'("../Mill_I/Hen15_plank/z0/out09.dat-0")')
        write(filefof,'("../Mill_I/Hen15_plank/z0/fof09.dat")')
        endif
        !----------------
        if(ncosmo.eq.4)then
        write(fileglx,'("guo_II/tablaglx.dat")')
        write(filerep,'("../Mill_II/Guo11_wmap1/z0/repetidas09.dat-0")')
        write(fileout,'("../Mill_II/Guo11_wmap1/z0/out09.dat-0")')
        write(filefof,'("../Mill_II/Guo11_wmap1/z0/fof09.dat")')
        endif
        !----------------
        if(ncosmo.eq.5)then
        !write(fileglx,'("guo_II_16/tablaglx.dat")')
        !write(filerep,'("guo_II_16/repetidas.dat-0")')
        !write(fileout,'("guo_II_16/out.dat-0")')
        !write(filefof,'("guo_II_16/fof.dat")')
        endif
        if(ncosmo.eq.7)then
        write(fileglx,'("hen_20/tablaglx.dat")')
        write(filerep,'("../Mill_I/Hen20_plank/z0/repetidas09.dat-0")')
        write(fileout,'("../Mill_I/Hen20_plank/z0/out09.dat-0")')
        write(filefof,'("../Mill_I/Hen20_plank/z0/fof09.dat")')
        endif
        !----------------
        if(ncosmo.eq.8)then
        write(fileglx,'("ayr_21/tablaglx.dat")')
        write(filerep,'("../Mill_I/Ayr21_plank/z0/repetidas09.dat-0")')
        write(fileout,'("../Mill_I/Ayr21_plank/z0/out09.dat-0")')
        write(filefof,'("../Mill_I/Ayr21_plank/z0/fof09.dat")')
        endif

        !----------------


c --------- abro archivos ----------------------------

        open(87,file=fileglx,status='unknown')

        open(88,file=filerep,status='old')
        rep=0
        do i=1,nmax
           read(88,*,end=10)col1,col2,t1,t2
           if(t1.le.t2)then
             rep(col2)=1
           else
             rep(col1)=1
           endif
        enddo
 10     nrep=i-1
        write(*,*)'repetidas',nrep
        close(88)
        
        open(89,file=filefof,status='old')
        do i=1,nmax
           read(89,*,end=11)ifila,igrut,inumt
           igru(ifila)=igrut
           inum(ifila)=inumt
        enddo
 11     nfof=i-1
        write(*,*)'filefof',nfof
        close(89)
       
        
        open(90,file=fileout,status='old')
        
        do i=1,nmax
                
                call percent(i,nmax,'escritura')
                if(ncosmo.eq.3.or.ncosmo.eq.7.or.ncosmo.eq.8)then
                   read(90,*,end=12)gID,ty,x,y,z,vx,vy,vz,stmass,bmass,
     &                              umag,gmag,rmag,imag,zmag,kmag,jmag
                else
                   read(90,*,end=12)gID,ty,x,y,z,vx,vy,vz,stmass,bmass,
     &                              umag,gmag,rmag,imag,zmag
                endif

                !print*,i,'leo'

                if(rep(i).EQ.1)go to 13

                dist=sqrt(x**2+y**2+z**2)

c                write(*,*)stmass,bmass
c                stop

                if(ncosmo.eq.3.or.ncosmo.eq.7.or.ncosmo.eq.8)then
                    write(87,78)gID,x,y,z,dist,vx,vy,vz,
     &                          umag,gmag,rmag,imag,zmag,kmag,jmag,
     &                          stmass,bmass,igru(i),inum(i)
                else
                    write(87,77)gID,x,y,z,dist,vx,vy,vz,
     &                          umag,gmag,rmag,imag,zmag,
     &                          stmass,bmass,igru(i),inum(i)
                endif


  

 13     end do

 12     nglx=i-1
        write(*,*)'galaxias',nglx
 77     format(i18,1x,7(f16.8,1x),5(f14.8,1x),2(e12.4,1x),2(i8,1x))
 78     format(i18,1x,7(f16.8,1x),7(f14.8,1x),2(e12.4,1x),2(i8,1x))
        close(87) 

        end
        include '../../../../subs/porcentajes.f'
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


