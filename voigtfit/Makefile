voigtfit.so: voigtfit.f95 evoigt.o mrqfit.o gaussj.o
	f2py -c -m voigtfit voigtfit.f95 evoigt.o mrqfit.o gaussj.o

evoigt.o: evoigt.f
	gfortran -c -O -fPIC evoigt.f

mrqfit.o: mrqfit.f
	gfortran -c -O -fPIC mrqfit.f

gaussj.o: gaussj.f
	gfortran -c -O -fPIC gaussj.f

clean: 
	rm mrqfit.o gaussj.o evoigt.o
