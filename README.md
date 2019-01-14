# TGE
This TGE package can be used to measure the angular power spectrum of the sky signal from radio interferometric visibilities. Details algorithm can be found in "https://arxiv.org/abs/1609.01732"

Required C packages:
(1) fftw3
(2) NRCP
(3) cfitsio
(4) gsl

Installation:
go to srcbin directory, run
make

Inputs:
#BCHAN
#ECHAN
#Umax(lamda)
#Umin(lamda)
#Ant_Dia(in m) 
#Theta0=0.6FWHM(arcmin) 
#f(taper)< 1 
#seed
#Nrea-Number of realizations
#flagLL (1=read LL 0=not read(-1.e8))
#fac (1=freq.scaling of the baselines, 0=no freq. scaling)
#Nbin(Number of bins)
#cutoff(some small number, 1.e-2 used)

