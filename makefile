LINKLIB=  -L/home/samir/astro/lib  -lfftw3 -lnrcp -lcfitsio -lm -lgsl -lgslcblas
INCLUDE=-I/home/samir/astro/include/

all: v2psfitsc binv2psfitscm


v2psfitsc: v2psfitsc.c read_fits_func.c v2psfitscfuncs.c beam.c
	gcc -o v2psfitsc $(INCLUDE) read_fits_func.c v2psfitscfuncs.c beam.c v2psfitsc.c $(LINKLIB)
	rm -rf *~

binv2psfitscm: binv2psfitscm.c
	gcc -o binv2psfitscm $(INCLUDE) binv2psfitscm.c $(LINKLIB)
	rm -rf *~
clean:
	rm -rf v2psfitsc binv2psfitscm
