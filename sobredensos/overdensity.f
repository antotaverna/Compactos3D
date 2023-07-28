!gfortran -o overdensity overdensity.f
!input: cosmologia 1, 2, 3, 5 ,7 o 8
!./overdensity input      
        program overdensity
        implicit none
        integer i,ii,n_snap,l
        integer is,nobj,nrep,lastsnap
        real  omega_m,over
        real  pi,delta_z0,denom,b0
        real  a,zz,b_z
        real,allocatable,dimension(:) :: z

        real evol
        external evol

c        integer*8 gID
c        real x,y,z,vx,vy,vz,umag,gmag,rmag,imag,zmag,stmass,bmass,sfr
c        real age,ubmag,gbmag,rbmag,ibmag,zbmag,dist

        character*80 string,filein,filein1,filein2
        character*80 fileout
        integer  iargc,ncosmo
        external iargc

        if (iargc().ne.1) then
         print *,'usage:  cosmo'
         print *,'Wmap_1/guo11 1 '
         print *,'Wmap_7/guo13 2 '
         print *,'Planck/hen20 3 '         
         print *,'Planck/ayr21 4 '         
         print *,'tng300 5 '         
         stop
        endif
        call getarg(1,string)
        read(string,*)ncosmo

        if(ncosmo.eq.1)then
                over = 25000. !sobredensidad en z=0
                omega_m = 0.25
                n_snap = 28   !(de z0 a z27)
                lastsnap = 64 !(de 0 a 63)
                filein='Mill_I/Guo11_wmap1/snapshots_wmap1.dat'
                fileout='../../files/guo_11/salida.dat'
        else if(ncosmo.eq.2)then
                over = 20000. !sobredensidad en z=0
                omega_m = 0.272 
                n_snap = 28
                lastsnap = 62 !(de 0 a 61)
                filein='Mill_I/Guo13_wmap7/snapshots_wmap7.dat'
                fileout='../../files/guo_13/salida.dat'
        else if(ncosmo.eq.3)then
                over = 20000. !sobredensidad en z=0
                omega_m = 0.315
                n_snap = 28
                lastsnap = 59 !(de 0 a 58) 
                filein='Mill_I/Hen20_plank/snapshots_plank.dat'
                fileout='../../files/hen_20/salida.dat'
        else if(ncosmo.eq.4)then
                over = 20000. !sobredensidad en z=0
                omega_m = 0.315
                n_snap = 28
                lastsnap = 59 !(de 0 a 58) 
                filein='Mill_I/Ayr21_plank/snapshots_plank.dat'
                fileout='../../files/ayr_21/salida.dat'        
        else if(ncosmo.eq.5)then
                !over = 20000. !sobredensidad en z=0
                !omega_m = 0.315
                !n_snap = 28
                !lastsnap = 59 !(de 0 a 58) 
                !filein='TNG300/XXX.dat'
                !fileout='../../files/tng300/salida.dat'        
        endif

        open(87,file=filein,status='old')
        open(90,file=fileout,status='unknown')

        allocate(z(n_snap))

        !over = 70000. !sobredensidad en z=0
        b0 = (over/evol(0.,omega_m))**(-1./3.)
        write(*,*)b0
        
        do i=1,100
             read(87,*,err=45)is,a,zz
             ii = lastsnap-is
             if(ii.le.n_snap)z(ii) = zz
        enddo
 45     close(87)

        !---------------------------------------------------------
        do i=1,1 !n_snap
           l=i-1
           if(l.lt.10)then
              if(ncosmo.eq.1)then
      write(filein1,'("Mill_I/Guo11_wmap1/z",i1,"/out09.dat-",i1)')l,l
      write(filein2,'("Mill_I/Guo11_wmap1/z",i1,"/repetidas09.dat-",i1)')l,l
              endif
              if(ncosmo.eq.2)then
      write(filein1,'("Mill_I/Guo13_wmap7/z",i1,"/out09.dat-",i1)')l,l
      write(filein2,'("Mill_I/Guo13_wmap7/z",i1,"/repetidas09.dat-",i1)')l,l
              endif
              if(ncosmo.eq.3)then
      write(filein1,'("Mill_I/Hen20_plank/z",i1,"/out09.dat-",i1)')l,l
      write(filein2,'("Mill_I/Hen20_plank/z",i1,"/repetidas09.dat-",i1)')l,l
              endif
              if(ncosmo.eq.4)then
      write(filein1,'("Mill_I/Ayr21_plank/z",i1,"/out09.dat-",i1)')l,l
      write(filein2,'("Mill_I/Ayr21_plank/z",i1,"/repetidas09.dat-",i1)')l,l
              endif
         endif

         !leo nlines de out.dat y repetidas.dat
         open(98,file=filein1,status='old')
         call count_lines(98,nobj)
         if(nobj.gt.15000000)write(*,*)'OJO DIMEN DE LO QUE SIGUE'
         open(99,file=filein2,status='old')
         call count_lines(99,nrep)
         !read(88,*)nobj
         !read(89,*)nrep             

         b_z = b0*(evol(z(i),omega_m))**(-1./3.)
         write(90,*)z(i),nobj,b_z,b_z**(-3.),nrep
         write(*,*)z(i),nobj,b_z,b_z**(-3.),nrep
        enddo
        !---------------------------------------------------------
        
        !close(88)
        !close(89)
        !close(90)
        close(98)      
        close(99)      

        end
        !b_z y b0 --> Ec. 3 paper Z & Euge 2014 
        !http://www.aanda.org/articles/aa/pdf/2014/01/aa22793-13.pdf 
c-----------------------------------------------------------------
        real function evol(z,omega_m)
        implicit none
        real z,omega_m,pi,delta_z0
        real factor

        factor = (1./omega_m-1.)*(1+z)**(-3) 
        pi = 4.*atan(1.)
        delta_z0 =18.*pi**2*(1+0.399*factor**0.941) 
        evol = 0.24*delta_z0/178+0.68

        return
        end
c-------------------------------
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
c-------------------------------
      subroutine getarg(arg,argt)
c     get a command line argument using gfortran
      implicit none
c input
c     number of argument to get
       integer arg
c output
c      string
       character*(*) argt
c executable
      call get_command_argument(arg,argt)
      return
      end
c---------------------------------------------------------------
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

