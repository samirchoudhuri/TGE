# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <gsl/gsl_sf_bessel.h>

extern double D;

double Beam(double pos,double freq)
{
  double D=45.,d=1.,nubyc,c;
  double be,arg;
  c=3.*pow(10,8.);//m/s
  
  nubyc=freq/c;//m^-1
  
  D=D*nubyc;
  d=d*nubyc;
  arg=M_PI*pos*D;
  if(pos==0.)
    {
      be=1.;
    }
  else
    {
      be=(2.)*(gsl_sf_bessel_J1(M_PI*pos*D)/(M_PI*pos*D));

      be=be*be;
    }
 
  //be=exp(-pos*pos/(9.8e-3*9.8e-3));
  return(be);
}
