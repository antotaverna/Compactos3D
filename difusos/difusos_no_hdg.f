c------gfortran -o difusos_no_hdg difusos_no_hdg.f
       program aisl
       implicit none 
       integer i,k,h,nlines,nn
       integer ngru,ngal

       !grupos
       integer j1  !nombre del grupo
       integer nmie,n3,nvec
       real xcm,ycm,zcm
       real vxcm,vycm,vzcm
       real sig_3d,rvir_bi,mass_bi
       real dist_media,dist_mediana
       real rabs1,rabs_n,deltamag
       real rmag_nvec,xvec,yvec,zvec
       real dmin,dmax
       real zc_min,zc_max,delta_zc
       real zs_min,zs_max,delta_zs

       integer nmatch
       real n_dif,n_hdg
       real mark(5000000)

       integer j2  !nombre del grupo

       integer difusos,marca

c-----------------------------------------
       external  iargc
       integer iargc,snap,ncosmo
       character*80 string
       character*80 filein1,filein2,filein0,filein3,filein4,fileout,fileout1
       
       if (iargc().ne.1) then
       print *,'Wmap_1/guo11 1 '
       print *,'Wmap_7/guo13 2 '
       print *,'Planck/hen15 3 '        
       print *,'Wmap_1/guo11-MII 4 '
       print *,'Wmap_1/guo11-MII_16 5 '
       print *,'Planck/hen20 7 '        
       stop
       endif

       call getarg(1,string)
       read(string,*)ncosmo 

c------------ FILENAMES -------------------------------   
       if(ncosmo.eq.1)then
           write(filein0,'("guo_11/tablaglx.dat")')
           write(filein1,'("guo_11/tablagru.dat")')
           write(filein4,'("guo_11/dif_eq_hdg.dat")')
           write(fileout,'("guo_11/difusos_no_hdg.dat")')
           !write(filein4,'("guo_11/dif_eq_wiens.dat")')
           !write(fileout,'("guo_11/difusos_no_wiens.dat")')
       endif
       !----------------
       if(ncosmo.eq.2)then
           write(filein0,'("guo_13/tablaglx.dat")')
           write(filein1,'("guo_13/tablagru.dat")')
           write(filein4,'("guo_13/dif_eq_hdg.dat")')
           write(fileout,'("guo_13/difusos_no_hdg.dat")')
       endif
       !----------------
       if(ncosmo.eq.3)then
           write(filein0,'("hen_15/tablaglx.dat")')
           write(filein1,'("hen_15/tablagru.dat")')
           write(filein4,'("hen_15/dif_eq_hdg.dat")')
           write(fileout,'("hen_15/difusos_no_hdg.dat")')
       endif
       !----------------
       if(ncosmo.eq.4)then
           write(filein0,'("guo_II/tablaglx.dat")')
           write(filein1,'("guo_II/tablagru.dat")')
           write(filein4,'("guo_II/dif_eq_hdg.dat")')
           write(fileout,'("guo_II/difusos_no_hdg.dat")')
       endif
       !----------------
       if(ncosmo.eq.5)then
           write(filein0,'("guo_II_16/tablaglx.dat")')
           write(filein1,'("guo_II_16/tablagru.dat")')
           write(filein4,'("guo_II_16/dif_eq_hdg.dat")')
           write(fileout,'("guo_II_16/difusos_no_hdg.dat")')
       endif
       !----------------
       if(ncosmo.eq.7)then
           write(filein0,'("hen_20/tablaglx.dat")')
           write(filein1,'("hen_20/tablagru.dat")')
           write(filein4,'("hen_20/dif_eq_hdg.dat")')
           write(fileout,'("hen_20/difusos_no_hdg.dat")')
       endif
       !----------------
c-----------------------------------------
       open(32,file=filein0,status='old')      
       open(30,file=filein1,status='old')
       open(34,file=filein4,status='old')      
       open(40,file=fileout,status='unknown')


       call count_lines(30,nlines)
       ngru = nlines
       print*,'numero de difusos ',ngru


       call count_lines(34,nlines)
       nmatch=nlines
       print*,'numero de difusos = HDG ',nmatch
       mark=0
       do i=1,nmatch
          read(34,*)n_dif,n_hdg
          mark(n_dif)=1
       enddo

       difusos=0
       do i=1,ngru
           read(30,*)j1,xcm,ycm,zcm,vxcm,vycm,vzcm,
     &                      sig_3d,nmie,rvir_bi,mass_bi,
     &                      rabs1,rabs_n,deltamag,n3



         !if(mark(j1).eq.1)goto 14 !elimino los HDG=Dif

         if(mark(j1).eq.0)then !me quedo con Difusos != HDG
            difusos=difusos+1
            write(40,77)j1,xcm,ycm,zcm,vxcm,vycm,vzcm,
     &                      sig_3d,nmie,rvir_bi,mass_bi,
     &                      rabs1,rabs_n,deltamag,n3,i

         endif       
 14    enddo  
       print*,'numero de grupos difusos no HDG', difusos, real(difusos)*100./real(ngru) 

       !IMPORTANTE
       !Las ultimas propiedades de los grupos aislados solo estan
       !bin calculadas para 1box.Cuando me paso a 8 boxes ya no me
       !sirven, porque cambian las coordenadas de los grupos.


 77     format(i8,1x,3(f10.5,1x),3(f12.5,1x),f10.3,1x,i8
     &         ,1x,f7.4,1x,e12.4,1x,3(f9.5,1x),i8,1x,i8)
 78     format(i8,15(1x,f17.8))      

       close(30)
       close(31)
       close(40)
       end program

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

c-----------------------------------------------------------------------
