# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include "read_fits_func.h"
# include "v2psfitscfuncs.h"

char INFITS[128],OUTFITS[128],input[128],OUTFITS_MG[128],OUTFITS_K2GG[128],OUTFITS_MG_K2GG[128];

extern int nrea;


int main(int argc, char *argv[])
{ 
  if(argc!=8)
    {
      printf("Usage: %s <input FITS file> <ouput FITS file> <input> <ouput FITS Mg file> <output data k2gg> <output Mg k2gg><FLAG (Note:0=Data Only; 1=MG Only; 2=Data+MG)>\n", argv[0]);
      return 1;
    }

  sscanf(argv[1],"%s",INFITS);  //Input: Input vis. fits file
  sscanf(argv[2],"%s",OUTFITS); //Input: Name of output fits file
  sscanf(argv[3],"%s",input);   //Input: Inputs
  sscanf(argv[4],"%s",OUTFITS_MG);   //Input: Name of output MG fits file
  sscanf(argv[5],"%s",OUTFITS_K2GG); 
  sscanf(argv[6],"%s",OUTFITS_MG_K2GG); 

  int mode;
  mode=atoi(argv[7]);
  printf("Running in mode %d.\n",mode);

  printf("\nReading from %s file\n",input);
  read_inputs(input);
  printf("Done reading inputs\n");
  
  //Check whether input fits file exists
  if(access(INFITS, F_OK)!=0){
    printf("Input File %s does not exists\n",INFITS);
    exit(0);
  }

  initialize();//Memory allocations and initialising seed for UAPS

  int flag=0;
  
  if(mode==0||mode==2)
    {
      printf("Starting correlation for data\n");
      correlate(flag);//Vis correlation for data
      printf("Data correlation done\n");
      flag=1;
      write_corr_fits(OUTFITS,OUTFITS_K2GG,flag);//Write vis. correlation for data
      printf("\nDone writing data\n");
    }
  
  if(mode==1||mode==2)
    {
      int ii;

      printf("\nStarting generating Mg\n");
      
      for(ii=0;ii<nrea;ii++)
	{
          Fill_UAPS();//Generate vis. for UAPS and fill them in *RRre,*RRim,*LLre,*LLim
	  correlate(ii);//Vis correlate for UAPS, add the values for nrea realizations
        }
      
      write_corr_fits(OUTFITS_MG,OUTFITS_MG_K2GG,nrea);//Write vis. correlation for MG
      printf("Done Generating Mg\n");

    }
  free_func();
  printf("\nThe End.\n");
}

