!gfortran -o readgroups2 readgroups2.f
!input1: cosmologia 1, 2, 3, 5, 7 o 8
!input2: redshift 0-27
!./readgroups2 input1 input2
!
!El archivo out.dat de Hen_15 y Hen_20 lee 2 lineas mas (kmag,jmag)
       program readdatos

        implicit none
        integer i,nrep,col1,col2,nmax,nglx
        parameter(nmax=15400000)
        integer rep(nmax),ifila,igrut,inumt,nfof
        integer igru(nmax),inum(nmax),ty
        integer ifila2,igrut2,inumt2,igru_d(nmax),inum_d(nmax)

        integer*8 gID
        real x,y,z,vx,vy,vz,umag,gmag,rmag,imag,zmag,kmag,jmag
        real stmass,bmass,sfr
        real age,ubmag,gbmag,rbmag,ibmag,zbmag,dist

        integer snap,ncosmo
        integer t1,t2

        external  iargc
        integer iargc
        character*80 string
        
        character*80 fileglx,filerep,fileout,filefof,filefof_d
c   
        integer ipcall,iplast,iphead,jc
        common /kpercent/ ipcall,iplast,iphead,jc

        if (iargc().ne.2) then
        print *,'usage:  ncosmo snap'
        print *,'Wmap_1/guo11 1 '
        print *,'Wmap_7/guo13 2 '
        print *,'Plank/hen15 3 '        
        print *,'Wmap_1/guo11-MII 4 '
        print *,'Wmap_1/guo11-MII_16 5 '
        print *,'Plank/hen20 7 '        
        print *,'Plank/ayr21 8 '        
        print *,'snap: z en el que estas parada (ejemplo: 5)'
        stop
        endif

        call getarg(1,string)
        read(string,*)ncosmo 
        call getarg(2,string)
        read(string,*)snap 

c------------ FILENAMES -------------------------------        
        if(ncosmo.eq.1)then
          if(snap.lt.10)then
        write(fileglx,'("guo_11/z",i1,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Guo11_wmap1/z",i1,"/repetidas09.dat-",i1)')snap,snap
      write(fileout,'("Mill_I/Guo11_wmap1/z",i1,"/out09.dat-",i1)')snap,snap
        write(filefof_d,'("Mill_I/Guo11_wmap1/z",i1,"/fof09.dat")')snap  !difusos
        write(filefof,'("guo_11/z",i1,"/fof09.dat")')snap  !denso
          else
        write(fileglx,'("guo_11/z",i2,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Guo11_wmap1/z",i2,"/repetidas09.dat-",i2)')snap,snap
      write(fileout,'("Mill_I/Guo11_wmap1/z",i2,"/out09.dat-",i2)')snap,snap
        write(filefof_d,'("Mill_I/Guo11_wmap1/z",i2,"/fof09.dat")')snap  !difusos
        write(filefof,'("guo_11/z",i2,"/fof09.dat")')snap
          endif
        endif
        !----------------
        if(ncosmo.eq.2)then
          if(snap.lt.10)then
        write(fileglx,'("guo_13/z",i1,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Guo13_wmap7/z",i1,"/repetidas09.dat-",i1)')snap,snap
      write(fileout,'("Mill_I/Guo13_wmap7/z",i1,"/out09.dat-",i1)')snap,snap
        write(filefof_d,'("Mill_I/Guo13_wmap7/z",i1,"/fof09.dat")')snap  !difusos
        write(filefof,'("guo_13/z",i1,"/fof09.dat")')snap
          else
        write(fileglx,'("guo_13/z",i2,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Guo13_wmap7/z",i2,"/repetidas09.dat-",i2)')snap,snap
      write(fileout,'("Mill_I/Guo13_wmap7/z",i2,"/out09.dat-",i2)')snap,snap
        write(filefof,'("guo_13/z",i2,"/fof09.dat")')snap
          endif
        endif
        !----------------
        if(ncosmo.eq.3)then
          if(snap.lt.10)then
        write(fileglx,'("hen_15/z",i1,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Hen15_plank/z",i1,"/repetidas09.dat-",i1)')snap,snap
      write(fileout,'("Mill_I/Hen15_plank/z",i1,"/out09.dat-",i1)')snap,snap
        write(filefof_d,'("Mill_I/Hen15_plank/z",i1,"/fof09.dat")')snap  !difusos
        write(filefof,'("hen_15/z",i1,"/fof09.dat")')snap
          else
        write(fileglx,'("hen_15/z",i2,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Hen15_plank/z",i2,"/repetidas09.dat-",i2)')snap,snap
      write(fileout,'("Mill_I/Hen15_plank/z",i2,"/out09.dat-",i2)')snap,snap
        write(filefof,'("hen_15/z",i2,"/fof09.dat")')snap
          endif
        endif
        !----------------
        if(ncosmo.eq.4)then
          if(snap.lt.10)then
        write(fileglx,'("guo_II/z",i1,"/tablaglx.dat")')snap
        write(filerep,'("Mill_II/Guo11_wmap1
     &/z",i1,"/repetidas09.dat-",i1)')snap,snap
        write(fileout,'("Mill_II/Guo11_wmap1
     &/z",i1,"/out09.dat-",i1)')snap,snap
        write(filefof_d,'("Mill_II/Guo11_wmap1/z",i1,"/fof09.dat")')snap  !difusos
        write(filefof,'("guo_II/z",i1,"/fof09.dat")')snap  !denso
          else
        write(fileglx,'("guo_II/z",i2,"/tablaglx.dat")')snap
        write(filerep,'("Mill_II/Guo11_wmap1
     &/z",i2,"/repetidas09.dat-",i2)')snap,snap
        write(fileout,'("Mill_II/Guo11_wmap1
     &/z",i2,"/out09.dat-",i2)')snap,snap
        write(filefof_d,'("Mill_II/Guo11_wmap1/z",i2,"/fof09.dat")')snap  !difusos
        write(filefof,'("guo_II/z",i2,"/fof09.dat")')snap
          endif
        endif
        !----------------
        if(ncosmo.eq.5)then
          if(snap.lt.10)then
        write(fileglx,'("guo_II/z",i1,"_16/tablaglx.dat")')snap
       write(filerep,'("guo_II/z",i1,"_16/repetidas09.dat-",i1)')snap,snap
        write(fileout,'("guo_II/z",i1,"_16/out09.dat-",i1)')snap,snap
        write(filefof,'("guo_II/z",i1,"_16/fof09.dat")')snap
          else
        write(fileglx,'("guo_II/z",i2,"_16/tablaglx.dat")')snap
       write(filerep,'("guo_II/z",i2,"_16/repetidas09.dat-",i2)')snap,snap
        write(fileout,'("guo_II/z",i2,"_16/out09.dat-",i2)')snap,snap
        write(filefof,'("guo_II/z",i2,"_16/fof09.dat")')snap
          endif
        endif
        !----------------
        if(ncosmo.eq.7)then
          if(snap.lt.10)then
        write(fileglx,'("hen_20/z",i1,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Hen20_plank/z",i1,"/repetidas09.dat-",i1)')snap,snap
        write(fileout,'("Mill_I/Hen20_plank/z",i1,"/out09.dat-",i1)')snap,snap
        write(filefof_d,'("Mill_I/Hen20_plank/z",i1,"/fof09.dat")')snap  !difusos
        write(filefof,'("hen_20/z",i1,"/fof09.dat")')snap
          else
        write(fileglx,'("hen_20/z",i2,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Hen20_plank/z",i2,"/repetidas09.dat-",i2)')snap,snap
        write(fileout,'("Mill_I/Hen20_plank/z",i2,"/out09.dat-",i2)')snap,snap
        write(filefof,'("hen_20/z",i2,"/fof09.dat")')snap
          endif
        endif
        !----------------
        if(ncosmo.eq.8)then
          if(snap.lt.10)then
        write(fileglx,'("ayr_21/z",i1,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Ayr21_plank/z",i1,"/repetidas09.dat-",i1)')snap,snap
        write(fileout,'("Mill_I/Ayr21_plank/z",i1,"/out09.dat-",i1)')snap,snap
        write(filefof_d,'("Mill_I/Ayr21_plank/z",i1,"/fof09.dat")')snap  !difusos
        write(filefof,'("ayr_21/z",i1,"/fof09.dat")')snap
          else
        write(fileglx,'("ayr_21/z",i2,"/tablaglx.dat")')snap
      write(filerep,'("Mill_I/Ayr21_plank/z",i2,"/repetidas09.dat-",i2)')snap,snap
        write(fileout,'("Mill_I/Ayr21_plank/z",i2,"/out09.dat-",i2)')snap,snap
        write(filefof,'("ayr_21/z",i2,"/fof09.dat")')snap
          endif
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
        
        open(89,file=filefof,status='old') !fof densos
        open(98,file=filefof_d,status='old') !fof difusos
        do i=1,nmax
           read(89,*,end=11)ifila,igrut,inumt
           read(98,*)ifila2,igrut2,inumt2
           igru(ifila)=igrut
           inum(ifila)=inumt
           igru_d(ifila)=igrut2
           inum_d(ifila)=inumt2
        enddo
 11     nfof=i-1
        write(*,*)'filefof',nfof
        close(89)
        close(98)
       
        
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
     &                          stmass,bmass,igru(i),inum(i),
     &                          igru_d(i),inum_d(i)
                else
                    write(87,77)gID,x,y,z,dist,vx,vy,vz,
     &                          umag,gmag,rmag,imag,zmag,
     &                          stmass,bmass,igru(i),inum(i),
     &                          igru_d(i),inum_d(i)
                endif


  

 13     end do

 12     nglx=i-1
        write(*,*)'galaxias',nglx
 77     format(i18,1x,7(f16.8,1x),5(f14.8,1x),2(e12.4,1x),4(i8,1x))
 78     format(i18,1x,7(f16.8,1x),7(f14.8,1x),2(e12.4,1x),4(i8,1x))
        close(87) 

        end
        include '/big4/users/subs/porcentajes.f'
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


