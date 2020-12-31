# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <gsl/gsl_sf_bessel.h>

extern double D;

double Beam(double pos,double freq)
{
  double nubyc,c,be;
  c=3.*pow(10,8.);//m/s

  nubyc=freq/c;//m^-1
  

  if(pos==0.)
    {
      be=1.;
    }
  else
    {
      be=(2.)*(gsl_sf_bessel_J1(M_PI*pos*D*nubyc)/(M_PI*pos*D*nubyc));

      be=be*be;
    }
 
  //be=exp(-pos*pos/(9.8e-3*9.8e-3));
  return(be);
}
