# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include <nr.h>
# include <nrutil.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include "read_fits_func.h"
# include "v2psfitscfuncs.h"

double Beam(double theta,double freq);
extern float *uu,*vv,**RRre,**RRim,**LLre,**LLim,del_chan,nu_chan0,chan0;
extern long nstokes,nchan,ncmplx,gcount;
extern char INFITS[128],OUTFITS[128],input[128],OUTFITS_MG[128];

double *GV,*GVS,*GVCre,*GVCim,GVCSQ=0.,*k2gg;
double *corrfact;


double theta_w,theta_0,D,theta_eff,f;
double Umax,Umin;
int chan1,chan2;
int nrea,flagLL;
double fac;

//For random number generation
unsigned long int seed;
gsl_rng *r;
double sigma=1.;
//done

double DD,dU;
long Nvis;
int n_ave,Ng,Nu,Nv,nbasln;
unsigned long *uindex,*vindex,total;

double L;//Resolution for UAPS image
int N;//No. of grid points for UAPS image


void printerror(int status)
{
  if (status)
    {
    fits_report_error(stderr, status);
    exit( status ); 
    }
}

void read_inputs(char *in_file)
{
  FILE *fp;
  fp=fopen(in_file,"r"); 
  fscanf(fp,"%d%d",&chan1,&chan2);
  fscanf(fp,"%lf%lf%lf%lf%lf%ld%d",&Umax,&Umin,&D,&theta_0,&f,&seed,&nrea);
  fscanf(fp,"%d%lf",&flagLL,&fac);
  fclose(fp);
}

void initialize()
{
  int ii;

  total=0;
  theta_0=(M_PI*theta_0)/(180.*60.);//theta_0 in rad
  theta_w=f*theta_0; // theta_w = f * theta_0
  theta_eff=f*theta_0/sqrt(1+f*f);

  printf("\ntheta_0=%e theta_w=%e theta_eff=%e\n",theta_0,theta_w,theta_eff);

  chan1=chan1-1; chan2=chan2-1;
  n_ave=chan2-chan1+1;  //no. of channels
  
  dU=0.265/(2.*theta_eff);// grid spacing
  DD=12.*dU;// uv separation upto which correlation is done={(ln 2)^(0.5)}/(2.*pi*theta_eff)

  printf("\n\nchan1-1=%d\tchan2-1=%d\tn_ave=%d\nUmax=%lf\tUmin=%lf\n\n",chan1,chan2,n_ave,Umax,Umin);
  printf("seed=%ld\tnrea=%d D=%lf\n\n",seed,nrea,D);
  printf("flagLL=%d\tfac=%f\n\n",flagLL,fac);
  
  printf("DD=%e(12dU)\tdU=%e\n",DD,dU);

  // read header, visibility data 
  //It will fill RRre[][],RRim[][],LLre[][],LLim[][], returns nbasln= no. of unflagged baselines

  nbasln=readfits(INFITS,Umax,Umin,chan1,chan2);//See read_fits_funcmfs.c
 
  printf("\nchan0=%e nu_chan0=%e del_chan=%e\n",chan0,nu_chan0,del_chan);
  printf("nbasln=%d nchan=%ld\n",nbasln,nchan);

  corrfact = (double*) calloc(nchan, sizeof(double));
  for(ii=0; ii<nchan; ii++)
    {    
      corrfact[ii] = 1.+ fac*(del_chan/nu_chan0)*(ii+1.+0.5-chan0);//lamda[chan0]/lambda[ii]
    }

  // dimensions  of array to store visibility correlation

  Ng =(int) ceil(((2.*Umax)-dU)/(2.*dU));
  Nu = 2*Ng +1;
  Nv = Ng + 1;
 
  //Nvis = Nu*Nv*n_ave;
  Nvis = Nu*Nv;

  printf("Nu=%d Nv=%d Nvis=%ld\n",Nu,Nv,Nvis);

  //Dimensions for UAPS image

  L=1./(2.*Umax);//L in radian
  N=(int)20*theta_0/L;//No. of grid points for UAPS
  if(N%2!=0) N=N+1;

  printf("L=%lf\tN=%d\n",L,N);

  /* index data according to u and v */
  total=(unsigned long)nbasln;
  uindex=(unsigned long*)calloc(nbasln,sizeof(unsigned long));
  vindex=(unsigned long*)calloc(nbasln,sizeof(unsigned long));

  //Index u,v values from uu[] & vv[] & store the indexed values in uindex[] & vindex[] in increasing order
  indexx(total,uu-1,uindex-1);
  indexx(total,vv-1,vindex-1);

  //for random number generators
  r= gsl_rng_alloc(gsl_rng_cmrg);
  gsl_rng_set (r, seed);
  //done
  printf("Initialization done\n");
}


void correlate(int fl)
{
  int mval;
  int ii,jj,kk,a,chan;
  double u1,u2,v1,v2;
  double diff1;
  int lmin,lmax;
  double wt1;
  double Re,Im;
  long index;

  if(fl==0)
    {   
      // arrays to store visibility correlation 
      GV = (double*) calloc(Nvis, sizeof(double));
      GVS= (double*) calloc(Nvis, sizeof(double));
      GVCre= (double*) calloc(Nvis, sizeof(double));
      GVCim= (double*) calloc(Nvis, sizeof(double));
      k2gg= (double*) calloc(Nvis, sizeof(double));
    }

  for(ii=0;ii<Nu;++ii)
    for(jj=0;jj<Nv;++jj)
      {
	index=jj*Nu+ii;
	GVCre[index]=0.;//Vc_re
	GVCim[index]=0.;//Vc_im
	GVS[index]=0.;//corr_re
      }

  lmin=0; lmax=0; 
  for(ii=0;ii<Nu;++ii) // set i grid coordinate 
    {
      // find left (lmin)  and right (lmax) limits  for baseline loop 
      while((lmin<nbasln)&&((ii-Ng)*dU-(uu[uindex[lmin]-1])>DD))	
	++lmin;
      lmax=lmin;
      while((lmax<nbasln)&&((uu[uindex[lmax]-1])-(ii-Ng)*dU<=DD))
	++lmax;
      // done select lmin and lmax 
      
      for(jj=0;jj<Nv;++jj) // set j grid coordinate 
	{
	  if(((ii-Ng)*(ii-Ng)+jj*jj)<=pow(Umax/dU,2.))
	    {
	      index=jj*Nu+ii;
	      for(kk=lmin;kk<lmax;++kk) // baseline loop  
		{
		  // identify baselines in relevant range and store index in mval
		  mval=uindex[kk]-1;
		  
		  u1=uu[mval];
		  v1=vv[mval];
		  
		  diff1=pow((ii-Ng)*dU-u1,2.)+pow(jj*dU-v1,2.);
		  
		  if(diff1 <= DD*DD)
		    {
		      for(chan=0;chan<n_ave;++chan)
			{
			  u1=uu[mval]*corrfact[chan+chan1];//modify u1,v1,diff1  with freq. correction
			  v1=vv[mval]*corrfact[chan+chan1];//here chan1=chan1-1
			  diff1=pow((ii-Ng)*dU-u1,2.)+pow(jj*dU-v1,2.);
			 
			  wt1=winf(diff1);
			  
			  if(RRre[mval][chan]>-1.e7)
			    {
			      Re=RRre[mval][chan];
			      Im=-1.*RRim[mval][chan];//as vis calculated with +1
			      
			      GVCre[index] +=(wt1*Re);//Vc_re
			      GVCim[index] +=(wt1*Im);//Vc_im
			      GVS[index] += (wt1*wt1*(Re*Re+Im*Im));//corr_re
			      k2gg[index]+=(wt1*wt1);
			    }
			  if(LLre[mval][chan]>-1.e7)
			    {
			      Re=LLre[mval][chan];
			      Im=-1.*LLim[mval][chan];//as vis calculated with +1
			      
			      GVCre[index] +=(wt1*Re);//Vc_re
			      GVCim[index] +=(wt1*Im);//Vc_im
			      GVS[index] += (wt1*wt1*(Re*Re+Im*Im));//corr_re
			      k2gg[index]+=(wt1*wt1);
			    }
			}//channel loop 
		    }// if
		}// end 'kk' loop
	    }//if loop
	}// for jj
    }// for ii
  
  printf("correlation done\n");
  
  for(ii=0;ii<Nu;++ii)
    for(jj=0;jj<Nv;++jj)
      {
	index=jj*Nu+ii;
	GVCSQ = (GVCre[index]*GVCre[index]+GVCim[index]*GVCim[index]);
	GV[index] += (GVCSQ - GVS[index]);//P(U) at grid point 
      }
}

double winf(double udif)
{
  double y;
  y=M_PI*pow(theta_w,2.)*exp(-1.*M_PI*M_PI*pow(theta_w,2.)*udif);
  return(y);
}

void write_corr_fits(char *OUT, char *OUTk2gg, int fl)
{
  fitsfile *fptr;
  int status=0,naxis=3;;
  char *CTYPE[3] = {"Nu", "Nv", "CHAN"};
  char keynam[6];
  double CRVAL[3],CDELT[3],CRPIX[3];
  long naxes[3];
  
  int ii;

  naxes[0]=(long)Nu;naxes[1]=(long)Nv;naxes[2]=1;
  CRVAL[0]=0.;CRVAL[1]=0.;CRVAL[2]=nu_chan0;
  CDELT[0]=dU;CDELT[1]=dU;CDELT[2]=del_chan;
  CRPIX[0]=(double) Ng;CRPIX[1]=1.;CRPIX[2]=chan0;
  
  fits_create_file(&fptr, OUT, &status);
  fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
  
  for(ii=0; ii<naxis; ii++)
    {
      if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_str(fptr, keynam, CTYPE[ii], " axis type", &status))
	printerror(status);
      
      if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_dbl(fptr, keynam, CRVAL[ii], 10, " ", &status))
	printerror(status);
      
      if(fits_make_keyn("CDELT", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_dbl(fptr, keynam, CDELT[ii],  9, " ", &status))
	printerror(status);
      
      if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_dbl(fptr, keynam, CRPIX[ii],  9, " ", &status))
	printerror(status);
    }

  for(ii=0;ii<Nvis;ii++)
   GV[ii]/=(double)fl;
 
  fits_write_img(fptr, TDOUBLE, 1, Nvis, GV, &status);
  fits_close_file(fptr, &status);

  naxes[0]=(long)Nu;naxes[1]=(long)Nv;naxes[2]=1;
  CRVAL[0]=0.;CRVAL[1]=0.;CRVAL[2]=nu_chan0;
  CDELT[0]=dU;CDELT[1]=dU;CDELT[2]=del_chan;
  CRPIX[0]=(double) Ng;CRPIX[1]=1.;CRPIX[2]=chan0;
  

  fits_create_file(&fptr, OUTk2gg, &status);
  fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
  
  for(ii=0; ii<naxis; ii++)
    {
      if(fits_make_keyn("CTYPE", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_str(fptr, keynam, CTYPE[ii], " axis type", &status))
	printerror(status);
      
      if(fits_make_keyn("CRVAL", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_dbl(fptr, keynam, CRVAL[ii], 10, " ", &status))
	printerror(status);
      
      if(fits_make_keyn("CDELT", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_dbl(fptr, keynam, CDELT[ii],  9, " ", &status))
	printerror(status);
      
      if(fits_make_keyn("CRPIX", ii+1, keynam, &status))
	printerror(status);
      if(fits_write_key_dbl(fptr, keynam, CRPIX[ii],  9, " ", &status))
	printerror(status);
    }
 
  for(ii=0;ii<Nvis;ii++)
    k2gg[ii]/=(double)fl;
 
  fits_write_img(fptr, TDOUBLE, 1, Nvis, k2gg, &status);
  fits_close_file(fptr, &status);

  free(k2gg);
  free(GVCre);
  free(GVCim);
  free(GVS);
  free(GV);
}

double  P_I(double nu,double U)//freq. in Hz
{
  double pu;
  double k_B=1.38e3,c=3.0e8;
  //pu=pow(2.*k_B*nu*nu/(c*c),2.)*A_150*pow(nu0/nu,2.0*alpha)*pow(1000./(2.*M_PI*U),betav);
  //pu=A_150*pow(1000./(2.*pi*U),betav);
  
  /////////////////////////////////////////////////////////////

  pu=pow(2.*k_B*nu*nu/(c*c),2.); // unit angular power spectrum 

  /////////////////////////////////////////////////////////////

  return(pu);
}

void Fill_UAPS()
{
  int i,j,k,index,index1,xdim,ydim,ia;
  double u,amp,nu0,nu,length,fac1;
  fftw_plan p,p1;
  fftw_complex *in,*reim;
  double *out,*img;
  double thetax,thetay,theta,be,spindex=-2.,delspindex=0.;

  fac1=L/(sqrt(2.)*N);length=N*L;
  ydim=(N/2+1);
  xdim=N;
  nu0=nu_chan0;
  printf("length=%e rad\n",length);

  out=(double*)calloc ((N*(N+2)),sizeof(double));
  in=(fftw_complex*)&out[0];
  img=(double*)calloc(N*N,sizeof(double));
  reim=(fftw_complex*)calloc ((N*(N/2+1)*nchan),sizeof(fftw_complex));

  p= fftw_plan_dft_c2r_2d (N, N, in,out, FFTW_ESTIMATE);
  p1= fftw_plan_dft_r2c_2d (N, N, out,in, FFTW_ESTIMATE);

  //Filling Fourier Components  
  //along axis (j-0 and j=N/2)
  for(j=0;j<ydim;j=j+N/2)
    for(i=1;i<N/2;++i)
      {
	// along + x 
	u=sqrt(1.*(i*i+j*j))/length;
	amp=fac1*sqrt(P_I(nu0,u));
	index=i*ydim+j;
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	
	// along -x 
	index1=(N-i)*ydim+j;
	in[index1][0]=in[index][0];
	in[index1][1]=-in[index][1];
	
      }
  // upper half plane excluding x axis
  
  for(i=0;i<xdim;++i)
    for(j=1;j<N/2;++j)
      {
	ia= (i>N/2) ? (N-i) : i ;
	u=sqrt(1.*(ia*ia+j*j))/length;
	amp=fac1*sqrt(P_I(nu0,u));
	index=i*ydim+j;
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
      }
  
  //4 points remain 
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	if(i+j==0) 
	  {
	    in[0][0]=0.0;
	    in[0][1]=0.0;
	  }
	else
	  {
	    u=(N/2.)*sqrt(1.*(i*i+j*j))/length;
	    amp=fac1*sqrt(P_I(nu0,u));
	    index=i*(N/2)*ydim+j*(N/2);
	    
	    in[index][0]=pow(-1.,(i*N/2+j*N/2))*amp*gsl_ran_gaussian(r,sigma);
	    in[index][1]=0.0;
	  }
      }
  
  // finished filling Fourier components

  fftw_execute(p);

  for(i=0;i<N;++i)
    for(j=0;j<N;++j)
      {
	index=i*(N+2)+j;
	index1=j*N+i;
	img[index1]=out[index];
      }

  for(k=0;k<nchan;k++)
    {
      nu=nu0+(k+1.+0.5-chan0)*del_chan;
      for(i=0;i<N;++i)
	for(j=0;j<N;++j)
	  {
	    index=i*(N+2)+j;
	    index1=j*N+i;
	    thetax= (i-N/2)*L;
	    thetay= (j-N/2)*L;
	    theta=sqrt((thetax*thetax)+(thetay*thetay));
	    be=Beam(theta,nu);
	    out[index]=be*img[index1]*pow((nu0/nu),spindex);
	  }
      fftw_execute(p1);
      for(i=0;i<N;++i)
        for(j=0;j<=N/2;++j)
          {
            index1=k*(N/2+1)*N+i*(N/2+1)+j;
            index=i*(N/2+1)+j;
            reim[index1][0]=pow(-1.,i+j)*in[index][0];//Image centre at (N/2,N/2)
            reim[index1][1]=pow(-1.,i+j)*(-1.*in[index][1]);//r2c use -1, to make it +1 exponent multiply with -1
          }
    }

  int chan;
  long group;
  double uuc,vvc,del_U;
  int ii,jj,ii1;
  
  del_U=(double)1/length;
  for(group=0;group<=nbasln;group++)
    {
      for(chan=0;chan<n_ave;chan++)
	{
	  uuc=uu[group]*corrfact[chan];
	  vvc=vv[group]*corrfact[chan];
	  if(abs(uuc)<(N*del_U/2.) && abs(vvc)<(N*del_U/2.))
	    {
	      ii1 = (int)roundf(uuc/del_U);
	      ii=(ii1<0) ? N+ii1 : ii1 ;
	      jj = (int)roundf(vvc/del_U);
	      
	      index1=chan*(N/2+1)*N+ii*(N/2+1)+jj;

	      if(RRre[group][chan]>-1.e7) RRre[group][chan]=reim[index1][0];
	      if(LLre[group][chan]>-1.e7) LLre[group][chan]=reim[index1][0];
	      if(RRim[group][chan]>-1.e7) RRim[group][chan]=reim[index1][1];
	      if(LLim[group][chan]>-1.e7) LLim[group][chan]=reim[index1][1];
	    }
	}
    }

  fftw_destroy_plan(p);
  fftw_destroy_plan(p1);
  fftw_free(out);
  fftw_free(reim);
}

void free_func()
{
  gsl_rng_free(r);
  free(uindex);
  free(vindex);
}
