/* cc -lm -o idenanalit idenanalit.c */
/* ./idenanalit ncosmo redshift  */

/* RECORDAR QUE EN C LA LECTURA DE ARCHIVOS NECESITA QUE ESPECIFIQUEN TOOOOODAS LAS COLUMNAS QUE TIENE EL ARCHIVO */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define SKIP fread(&dummy, sizeof(dummy), 1, pf);

int grid(float *x,float *y,float *z, float rlbox, int ngrid, int nobj, 
         int *lirst, int *ll);

int busv(float *x, float *y, float *z, int ic, float rlbox, int ngrid, int nobj,
         int *lirst, int *ll, int *grupo, float r0, int iper, int *listvec, 
         int *nvec);

int iden(float *x, float *y, float *z, float rlbox, float r0, int nobj, 
         int iper, int *grupo);

int main(int argc, char **argv)
{
  int   i,nid,nnid,t1,t2;
  int   *id;
  int   *marca,*marca99;
  int   nobj,ncont,nrep,nmag;
  int   ii;
  char  idg[20];
  float *x,*y,*z;
  float rlbox;
  int iper;
  int *grupo;
  int *numero;
  double r0;
  char idglx[20];
  int idd;
  int type;
  float rx,ry,rz,jreal;
  float vx,vy,vz,umag,gmag,rmag,imag,zmag,kmag,jmag;
  float stmass,bmass,sfr,age,ubmag,gbmag,rbmag,ibmag,zbmag;

  char filename1[200],filename2[200],filename3[200],filename[200];
  FILE  *pf,*pf1;
  int dummy, n, m, k, pc;
  float r[3];
  int redshift,ncosmo;
  int nn,nr,*nsal,*nsal2;
  float zzt,bb,*bz,dd;
  int n_snap;
  float xmax,ymax,zmax; 
	
/******************** INPUT **********************************/
 if (argc == 3){
   ncosmo = atoi(argv[1]); // use command line parameter; // use default from ANAL_PARAMS.H
  }else{
   fprintf(stderr, "USAGE: idenanalit <redshift> <ncosmo>..\n");
   return -1;
 }

 if (argc == 3){
   redshift = atoi(argv[2]); // use command line parameter; // use default from ANAL_PARAMS.H
  }else{
   fprintf(stderr, "USAGE: idenanalit <redshift> ..\n");
   return -1;
 }

 if(ncosmo == 1)printf("cosmo:Guo11_wmap1=%d\n",ncosmo);
 if(ncosmo == 2)printf("cosmo:Guo13_wmap7=%d\n",ncosmo);
 if(ncosmo == 3)printf("cosmo:Hen15_plank=%d\n",ncosmo);
 if(ncosmo == 4)printf("cosmo:=Guo11_wmap1-MII=%d\n",ncosmo);
 if(ncosmo == 5)printf("cosmo:=Guo11_wmap1-MII_16=%d\n",ncosmo);
 if(ncosmo == 6)printf("cosmo:=Guo11_wmap1-MII_16difu=%d\n",ncosmo);
 if(ncosmo == 7)printf("cosmo:=Hen20_plank=%d\n",ncosmo);
 if(ncosmo == 8)printf("cosmo:=Ayr21_plank=%d\n",ncosmo);
 printf("z=%d\n",redshift);

/*************************************************************/

 if(ncosmo == 1){
      sprintf(filename1,"guo_11/salida.dat");
      sprintf(filename2, "Mill_I/Guo11_wmap1/z%i/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "Mill_I/Guo11_wmap1/z%i/out09.dat-%i", redshift, redshift);
      sprintf(filename, "guo_11/z%i/fof09.dat", redshift);
      rlbox = 500000;/*en kpc*/
 }
 if(ncosmo == 2){
      sprintf(filename1,"guo_13/salida.dat");
      sprintf(filename2, "Mill_I/Guo13_wmap7/z%i/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "Mill_I/Guo13_wmap7/z%i/out09.dat-%i", redshift, redshift);
      sprintf(filename, "guo_13/z%i/fof09.dat", redshift);
      rlbox = 500000;/*en kpc*/
 }
 if(ncosmo == 3){
      sprintf(filename1,"hen_15/salida.dat");
      sprintf(filename2, "Mill_I/Hen15_plank/z%i/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "Mill_I/Hen15_plank/z%i/out09.dat-%i", redshift, redshift);
      sprintf(filename, "hen_15/z%i/fof09.dat", redshift);
      rlbox = 480279;/*en kpc*/
 }
 if(ncosmo == 4){
      sprintf(filename1,"guo_II/salida.dat");
      sprintf(filename2, "Mill_II/Guo11_wmap1/z%i/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "Mill_II/Guo11_wmap1/z%i/out09.dat-%i", redshift, redshift);
      sprintf(filename, "guo_II/z%i/fof09.dat", redshift);
      rlbox = 100000;/*en kpc*/
 }
/* corregida para correr sobredensidad 5000 a z=0*/
 if(ncosmo == 5){
      sprintf(filename1,"guo_II/salida_16.dat");
      sprintf(filename2, "guo_II/z%i_16/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "guo_II/z%i_16/out09.dat-%i", redshift, redshift);
      sprintf(filename, "guo_II/z%i_16/fof09.dat", redshift);
      rlbox = 100000;/*en kpc*/
 }
 if(ncosmo == 6){
      sprintf(filename1,"guo_II/salida_16difu.dat");
      sprintf(filename2, "guo_II/z%i_16/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "guo_II/z%i_16/out09.dat-%i", redshift, redshift);
      sprintf(filename, "difusos/guo_II_16/fof09.dat");
      rlbox = 100000;/*en kpc*/
 }
if(ncosmo == 7){
      sprintf(filename1,"hen_20/salida.dat");
      sprintf(filename2, "Mill_I/Hen20_plank/z%i/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "Mill_I/Hen20_plank/z%i/out09.dat-%i", redshift, redshift);
      sprintf(filename, "hen_20/z%i/fof09.dat", redshift);
      rlbox = 480279;/*en kpc*/
 }
if(ncosmo == 8){
      sprintf(filename1,"ayr_21/salida.dat");
      sprintf(filename2, "Mill_I/Ayr21_plank/z%i/repetidas09.dat-%i", redshift,redshift);
      sprintf(filename3, "Mill_I/Ayr21_plank/z%i/out09.dat-%i", redshift, redshift);
      sprintf(filename, "ayr_21/z%i/fof09.dat", redshift);
      rlbox = 480279;/*en kpc*/
 }



 if(!(pf1=fopen(filename1,"r")))
  {
    printf("can't open file `%s`\n",filename1);
    exit(0);
  }
  printf("open file `%s`\n",filename1);

  n_snap = 28;
  nsal         = malloc(n_snap*sizeof(int));
  nsal2        = malloc(n_snap*sizeof(int));
  bz           = malloc(n_snap*sizeof(float));

  for(n=0;n<n_snap;n++)
  { fscanf(pf1, "%f %d %f %f %d\n",&zzt,&nn,&bb,&dd,&nr);
      nsal[n]=nn;
      bz[n]=bb;
      nsal2[n]=nr;
      //printf("%d %d %f\n",n,nn,bb);
  }
  fclose(pf1);

  nobj=nsal[redshift];
      printf("n_out %d\n",nobj);
  nrep=nsal2[redshift];
      printf("n_rep %d\n",nrep);

  if(!(pf1=fopen(filename2,"r")))
  {
    printf("can't open file `%s`\n",filename2);
    exit(0);
  }
  printf("open file `%s`\n",filename2);

  marca         = malloc(nobj*sizeof(int));
  for(n=0;n<nobj;n++)
     { marca[n]=0;
     }
    for(n=0;n<nrep;n++)
    { fscanf(pf1, "%d %d %d %d\n",&ii, &nid, &t1, &t2);
      if(t1 <= t2){
      	   nnid=nid-1;	    
      }else{
      	   nnid=ii-1;	    
      }
      marca[nnid]=1;
    }
  fclose(pf1);
  

  if(!(pf=fopen(filename3,"r")))
  {
    printf("can't open file `%s`\n",filename3);
    exit(0);
  }
  printf("paso el open3\n");

  id         = malloc(nobj*sizeof(int));
  x          = malloc(nobj*sizeof(float));
  y          = malloc(nobj*sizeof(float));
  z          = malloc(nobj*sizeof(float));
  grupo      = malloc(nobj*sizeof(int));
  numero     = malloc(nobj*sizeof(int));

  //rlbox = 500000;/*en kpc*/
  ncont=0;
      printf("%d\n",nobj); 

      xmax=-1.;
      ymax=-1.;
      zmax=-1.;

    for(n=0;n<nobj;n++)
    { 
      if(ncosmo == 1){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag);}
      if(ncosmo == 2){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag);}
       /*para hen_15 leo 2 col mas */
      if(ncosmo == 3){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag,&kmag,&jmag);}
      if(ncosmo == 4){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag);}
      if(ncosmo == 5){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag);}
      if(ncosmo == 6){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag);}
       /*para hen_20 leo 2 col mas */
      if(ncosmo == 7){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag,&kmag,&jmag);}
      if(ncosmo == 8){fscanf(pf, "%s %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",&idglx,&type,&rx,&ry,&rz,&vx,&vy,&vz,&stmass,&bmass,&umag,&gmag,&rmag,&imag,&zmag,&kmag,&jmag);}

      /*printf("%d %s\n",n,idglx); */
     jreal = (float)n/(float)1000000-(float)((int)((float)n/(float)1000000));
    /*if(jreal == 0.)
    {
      printf("%d\n",n);
    }*/
    /*  if(marca[idd] == 0 && marca99[idd] == 0){*/
      idd=n+1;
     if(marca[n] == 0){
      ncont=ncont+1;
      id[ncont]=idd;
      r[0]=rx ; 
      r[1]=ry ; 
      r[2]=rz ;
      r[0]*=1000. ;
      r[1]*=1000. ;
      r[2]*=1000. ;
      if(r[0] == rlbox ) r[0] = r[0]-1.;
      if(r[1] == rlbox ) r[1] = r[1]-1.;
      if(r[2] == rlbox ) r[2] = r[2]-1.;
      if(r[0] < 0.)    r[0] += rlbox ;
      if(r[0] > rlbox) r[0] -= rlbox ;
      if(r[1] < 0.)    r[1] += rlbox ;
      if(r[1] > rlbox) r[1] -= rlbox ;
      if(r[2] < 0.)    r[2] += rlbox ;
      if(r[2] > rlbox) r[2] -= rlbox ;
      if(xmax < r[0]) xmax = r[0];
      if(ymax < r[1]) ymax = r[1];
      if(zmax < r[2]) zmax = r[2];
      x[ncont]=r[0];
      y[ncont]=r[1];
      z[ncont]=r[2];
/*	printf("%f %f %f\n",x[ncont],y[ncont],z[ncont]);*/
      } 
    }

  fclose(pf);
  nobj=ncont;
      printf("%d\n",nobj); 
      printf("%f %f %f\n",xmax,ymax,zmax); 

  iper=1;/*periodico*/
  /* yose paper I: r0 = 0.132177*rlbox/pow((float)nobj,0.33333);  *//*0.14 para identificar halos de galaxias */

//  r0 = bz[redshift]*rlbox/pow((float)nobj,0.33333); 
  r0 = 90. ;    //0.14 para identificar halos de galaxias
  printf("nobj rlbox r0 %d %f %f %f\n",nobj,rlbox,r0,bz[redshift]);

  iden(x,y,z,rlbox,r0,nobj,iper,grupo);

  for(i=0;i<nobj;i++)
    numero[i]=0;
  for(i=0;i<nobj;i++)
    if(grupo[i]>0)numero[grupo[i]]++;

//  pf=fopen("fof.dat","w+");
  pf=fopen(filename,"w+");
  for( i = 1 ; i < nobj+1; i++ )
    /*fprintf(pf,"%d %d %f %f %f %f %f %f %d %d\n",id[i],idglx[i],x[i],y[i],z[i],vx[i],vy[i],vz[i],grupo[i],numero[grupo[i]]);*/
/*{ if( numero[grupo[i]] < 4 ){
       grupo[i]=0 ;
       numero[grupo[i]]=0 ;  
}*/
    fprintf(pf,"%d %d %d\n",id[i],grupo[i],numero[grupo[i]]);
  fclose(pf);
}

/*********************** FUNCIONES ********************************/

int iden(float *x, float *y, float *z, float rlbox, float r0, int nobj, 
         int iper, int *grupo)
{
/*se asume que las posisiones estan entre 0 y rlbox*/
  int   i, j;
  int   *lirst, *ll;
  int   *listvec, *listvecnew;
  int   ngrid;
  int   nvec, nvecvec;
  int   ngrupos;
  float jreal;

  ngrid=(int)(rlbox/r0);
  if(ngrid>256)ngrid=256;
  printf("r0 ngrid %g %d\n",r0,ngrid);
  lirst      = malloc(ngrid*ngrid*ngrid*sizeof(int));
  ll         = malloc(nobj*sizeof(int));
  listvec    = malloc(nobj*sizeof(int));
  listvecnew = malloc(nobj*sizeof(int));

  grid(x,y,z,rlbox,ngrid,nobj,lirst,ll);

/*empieza la identificacion*/
  for( i = 0 ; i < nobj ; i++ )
    grupo[i]=0;
  
  ngrupos = 0;
/*  pf=fopen("identest.dat","w+");*/
  for( i = 0 ; i < nobj; i++ )
  {
    jreal = (float)i/(float)1000000-(float)((int)((float)i/(float)1000000));
    if(jreal == 0.)
    {
      printf("%d\n",i);
    } 
    if( grupo[i] != 0 ) continue; /*salta a la siguiente*/
    //printf("%d de %d\n",i,nobj);
    nvec=0;
    busv(x,y,z,i,rlbox,ngrid,nobj,lirst,ll,grupo,r0,iper,listvec,&nvec);
    //printf("nvec %f %f %f %d\n",x[i],y[i],z[i],nvec);
    if(nvec!=0)
    {
      ngrupos++;
      grupo[i]=ngrupos;
    }
    for( j = 0 ; j < nvec ; j++ )
      grupo[listvec[j]]=ngrupos;
	

    do
    {
      nvecvec=0;
      for( j = 0 ; j < nvec ; j++ )
      {
        int k;
        busv(x,y,z,listvec[j],rlbox,ngrid,nobj,lirst,ll,grupo,r0,iper,
             listvecnew,&nvecvec);
        for(k=0;k<nvecvec;k++)
          grupo[listvecnew[k]]=ngrupos;
      }

      if( nvecvec == 0 ) break;
      for( j = 0 ; j < nvecvec ; j++ )
        listvec[j]=listvecnew[j];
      nvec=nvecvec;
/*      printf("      nvecvec %d\n",nvecvec);*/
/*      for(j=0;j<nvecvec;j++)
      {
        int ind;
        ind=listvecnew[j];
        fprintf(pf,"%f %f %f %d %d\n",x[ind],y[ind],z[ind],ind,grupo[ind]);
      }
      exit(0);*/
    } while( 1 ); /*fin del do de los vecinos de los vecinos*/
  } 
/*  fclose(pf);*/
/*termino la identificacion*/
}

/*********************************/
int grid(float *x,float *y,float *z, float rlbox, int ngrid, int nobj, 
         int *lirst, int *ll)
{
  int i, j, k;
  int ix, iy, iz;
  float fac;

  fac = (float)ngrid/rlbox ;

  printf("grid\n");
  for( i = 0 ; i < ngrid ; i++)
    for( j = 0 ; j < ngrid ; j++)
      for( k = 0 ; k < ngrid ; k++)
        lirst[(i * ngrid + j) * ngrid + k] = 0 ;

  for( i = 0 ; i < nobj ; i++ )  
    ll[i] = 0 ;

  for( i = 0 ; i < nobj ; i++ )
  {
    ix=(int)(x[i]*fac);
    iy=(int)(y[i]*fac);
    iz=(int)(z[i]*fac);
/*printf("%d %d.........%d\n",ix,iy,iz);
printf("%f %f.........%f\n",x[i],y[i],z[i]);*/
    lirst[(ix * ngrid + iy) * ngrid + iz] = i ;
  }  

  for( i = 0 ; i < nobj ; i++ )
  {
    ix=(int)(x[i]*fac);
    iy=(int)(y[i]*fac);
    iz=(int)(z[i]*fac);
    ll[lirst[(ix * ngrid + iy) * ngrid + iz]] = i ;
    lirst[(ix * ngrid + iy) * ngrid + iz] = i ;
  }  

  printf("fin grid\n");
}

/*********************************/
int busv(float *x, float *y, float *z, int ic, float rlbox, int ngrid, int nobj,
         int *lirst, int *ll, int *grupo, float r0, int iper, int *listvec, 
         int *nvec)
{
  int ixc, iyc, izc;
  int ixci, iyci, izci;
  int ixcf, iycf, izcf;
  int ix, iy, iz;
  int ixx, iyy, izz;
  int ibox;
  float xx, yy, zz, dis;
  int i;
  float rlbox2;

  rlbox2=rlbox/2.0;
  ixc=(int)(x[ic]/rlbox*ngrid);
  ixci=ixc-1;
  ixcf=ixc+1;
  iyc=(int)(y[ic]/rlbox*ngrid);
  iyci=iyc-1;
  iycf=iyc+1;
  izc=(int)(z[ic]/rlbox*ngrid);
  izci=izc-1;
  izcf=izc+1;
  if ( iper == 0 ) /*no periodico*/
  {
    if( ixci < 0 ) ixci = 0;
    if( iyci < 0 ) iyci = 0;
    if( izci < 0 ) izci = 0;
    if( ixcf >= ngrid ) ixcf = ngrid-1;
    if( iycf >= ngrid ) iycf = ngrid-1;
    if( izcf >= ngrid ) izcf = ngrid-1;
  } 
  else if ( iper != 1 )
  {
    printf("la periodicidad debe ser 0 o 1\n");
    exit(0);
  }
  for( ixx = ixci ; ixx <= ixcf ; ixx++)
  {
    ix=ixx;
    if(iper == 1)
    {
      if(ix >= ngrid) ix = ix - ngrid;
      if(ix < 0) ix = ix + ngrid;
    }
    for( iyy = iyci ; iyy <= iycf ; iyy++)
    {
      iy=iyy;
      if(iper == 1)
      {
        if(iy >= ngrid) iy = iy - ngrid;
        if(iy < 0) iy = iy + ngrid;
      }
      for( izz = izci ; izz <= izcf ; izz++)
      {
        iz=izz;
        if(iper == 1)
        {
          if(iz >= ngrid) iz = iz - ngrid;
          if(iz < 0) iz = iz + ngrid;
        }

        ibox=(ix * ngrid + iy) * ngrid + iz ;
        i=lirst[ibox];
        if(i==-1)continue;
        do
        {
	  /*if(ix == 203 && iy == 157 && iz == 208 && i != 0)
	  {
	  	printf("dis %d %d %d %f %f %f\n",ibox,ll[i],i,x[i],x[ic],rlbox);
	  }*/
          i=ll[i];
          xx=x[i]-x[ic];
          yy=y[i]-y[ic];
          zz=z[i]-z[ic];
          if( iper == 1 )
          {
            if( xx > rlbox2 ) xx = xx - rlbox;
            if( yy > rlbox2 ) yy = yy - rlbox;
            if( zz > rlbox2 ) zz = zz - rlbox;
            if( xx < -rlbox2 ) xx = xx + rlbox;
            if( yy < -rlbox2 ) yy = yy + rlbox;
            if( zz < -rlbox2 ) zz = zz + rlbox;
          }
          dis=sqrt(xx*xx+yy*yy+zz*zz);
          if( dis < r0 && grupo[i] == 0 && i != ic )
          {
            listvec[*nvec]=i;
            *nvec=*nvec+1;
	    //printf("muchos vecinos %d\n",*nvec);
          }
        } while( i != lirst[ibox] ); /*fin laso particulas del grid*/
      } /*fin izz*/
    } /*fin iyy*/
  } /*fin ixx*/
}
