# TGE
This Tapered Gridded Estimator (TGE) package can be used to estimate the angular power spectrum of the sky signal from radio interferometric visibilities. Details algorithm can be found in "https://arxiv.org/abs/1609.01732"

Required C packages:
(1) fftw3
(2) NRCP
(3) cfitsio
(4) gsl

Installation:
make

Inputs:
#BCHAN (start channel number)

#ECHAN (end channel number)

#Umax(lamda) (max baseline length in lambda unit)

#Umin(lamda)  (min baseline length in lambda unit)

#Ant_Dia(in m) (Antenna diameter for primary beam calculation in m)

#Theta0=0.6FWHM(arcmin)

#f(taper)< 1 (set tapering window)

#seed (for random realizations)

#Nrea-Number of realizations

#flagLL (1=read LL polarization, 0=not read(-1.e8))

#fac (1=freq.scaling of the baselines, 0=no freq. scaling)

#Nbin(Number of bins)

#cutoff(some small number, 1.e-2 used)

After giving the input parameters:

sh run.sh

