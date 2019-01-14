// This program reads the gridded visibilities,and make image write in FITS format 
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>

double CRVAL[3],CDELT[3],CRPIX[3];
long naxes[3];
int naxis=3,nfound;

int main(int argc, char *argv[])
{ 
  char INFITS[128],CONST[128];
  double Umax,Umin,kv,cutoff;
  int Nbin;
  FILE *fp;

  if(argc!=5){
    printf("Usage: %s <input FITS file> <input file> <outputfile> <const FITS file>\n", argv[0]);
    return 1;
  }
  
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[4],"%s",CONST);
  
  //reading input parameter
  fp=fopen(argv[2],"r");
  fscanf(fp,"%*d%*d");
  fscanf(fp,"%lf%lf%*f%*f%*f%*d%*d",&Umax,&Umin);
  fscanf(fp,"%*d%*f");
  fscanf(fp,"%d%lf",&Nbin,&cutoff);
  fclose(fp);
  
  printf("Umax=%.2f Umin=%.2f Nbin=%d  cutoff=%e\n",Umax,Umin,Nbin,cutoff);
  kv=(1.*Nbin)/log10(Umax/Umin);
  
  fitsfile *fptr,*fptr1;
  int Nu,Nv,Ng,Nuv,n_ave,chan;
  long Nvis;
  int status=0,anynull;
  long fpixel[3],lpixel[3],inc[3]={1,1,1};
  double nulval=0;  
  double *GV,*GVconst,nullval;
  double UvalGrid,*binval,*uval,*weight,*num,Ag2,Bg,Cg2,wg;
  int NUGrid,ii,jj,ii1,index;
  
  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }
  
  if(access(CONST, F_OK)!=0){
    printf("Input File %s does not exists\n",CONST);
    exit(0);
  }
   
  //read input header  
  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_read_keys_lng(fptr,"NAXIS",1,naxis,naxes,&nfound,&status);
  fits_read_keys_dbl(fptr,"CDELT",1,naxis,CDELT,&nfound,&status);
  fits_read_keys_dbl(fptr,"CRVAL",1,naxis,CRVAL,&nfound,&status);
  fits_read_keys_dbl(fptr,"CRPIX",1,naxis,CRPIX,&nfound,&status);
  fits_close_file(fptr,&status);

  printf("dU=%lf naxes vis={%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2]);
  Nu=(int) naxes[0]; Ng=(int) CRPIX[0]; Nv=(int) naxes[1]; n_ave=(int) naxes[2];
  Nuv= Nu*Nv;
  printf("Nu=%d Nv=%d Ng=%d Nuv=%d n_ave=%d\n",Nu,Nv,Ng,Nuv,n_ave);
  
  //read input header  
  fits_open_file(&fptr,CONST,READONLY,&status);
  fits_read_keys_lng(fptr,"NAXIS",1,naxis,naxes,&nfound,&status);
  fits_read_keys_dbl(fptr,"CDELT",1,naxis,CDELT,&nfound,&status);
  fits_read_keys_dbl(fptr,"CRVAL",1,naxis,CRVAL,&nfound,&status);
  fits_read_keys_dbl(fptr,"CRPIX",1,naxis,CRPIX,&nfound,&status);
  fits_close_file(fptr,&status);
    
  printf("dU=%lf naxes vis={%ld,%ld,%ld}\n",CDELT[0],naxes[0],naxes[1],naxes[2]); 
  Nu=(int) naxes[0];Nv=(int) naxes[1];
  printf("Nu=%d Nv=%d\n",Nu,Nv);

 
  GV = (double*) calloc(Nuv, sizeof(double));
  GVconst = (double*) calloc(Nuv, sizeof(double));
  binval=(double*)calloc(Nbin,sizeof(double));
  weight=(double*)calloc(Nbin,sizeof(double));
  uval=(double*)calloc(Nbin,sizeof(double));
  num=(double*)calloc(Nbin,sizeof(double));
  
  fits_open_file(&fptr,INFITS,READONLY,&status);
  fits_open_file(&fptr1,CONST,READONLY,&status);
  fp = fopen(argv[3],"w");  
  
  for(chan=0;chan<n_ave;++chan)
    { 
      printf("chan=%d\n",chan); 
      fpixel[0]=1;fpixel[1]=1;fpixel[2]=chan+1;
      lpixel[0]=(long)Nu;lpixel[1]=(long)Nv;lpixel[2]=chan+1;
      fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,&nulval,GV,&anynull,&status);
      fits_read_subset(fptr1,TDOUBLE,fpixel,lpixel,inc,&nulval,GVconst,&anynull,&status);
  
      for(ii=0;ii<Nbin;++ii)
      	{
      	  uval[ii]=0.;
      	  binval[ii]=0.;
      	  weight[ii]=0.;
	  num[ii]=0.;
      	}
      
      for(jj=0;jj<Nv;jj++)
      	for(ii=0;ii<Nu;ii++)
      	  {
      	    index=jj*Nu+ii;
      	    ii1=ii-Ng;
      	    UvalGrid=CDELT[0]*sqrt(ii1*ii1+jj*jj);
      	    if(UvalGrid> Umin && UvalGrid<= Umax)
      	      {
		NUGrid=(int)floor(kv*log10(UvalGrid/Umin));
      		NUGrid = (NUGrid>0) ? NUGrid:0;
      		Ag2=GV[index];
		Cg2=GVconst[index];
      		//wg=GVconst[index];
		wg=1.;
		if(Cg2>cutoff && NUGrid<Nbin)
		  {
		    binval[NUGrid]+=(wg*Ag2/Cg2);
		    uval[NUGrid]+=(UvalGrid*wg);
		    weight[NUGrid]+=wg;
		    num[NUGrid]+=1.;
		  }
	      }
      	  }
      
      for(ii=0;ii<Nbin;++ii)
      	{
      	  if(weight[ii]>0.)
      	    fprintf(fp,"%e %e %e\n",uval[ii]/weight[ii],binval[ii]/weight[ii],num[ii]);
	  else
      	    fprintf(fp,"%e %e %e\n",0.,-1.e10,0.);
      	}
      fprintf(fp,"\n");
    }
  
  fits_close_file(fptr,&status);
  fits_close_file(fptr1,&status);  
  fclose(fp);
}
