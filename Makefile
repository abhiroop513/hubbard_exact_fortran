F90S= gfortran

#FLAGS=  -O3 -fcheck=bounds
FLAGS=  -O3  -cpp -ffree-form -ffast-math -ftree-vectorize -funroll-loops -ffree-line-length-none 
MATHLIBS= -lblas -llapack
#MATHLIBS= -framework Accelerate -bind_at_load
FPPFLAGS=

gfortran -fno-range-check -c paramod.f90
gfortran -fno-range-check -c interface.f90
./f92 gendet-amlt
./f92 matpppinit
./f92 makeinit
./f6  # diagonly
./f7  # lanc_dsf (for Dynamical Structure Factor)
gfortran -fcheck=bounds -mcmodel=large diagonly.f90 -L/home/abhiroop/hubbard_package -llapack -L/home/abhiroop/hubbard_package -lrefblas -o diagonly.x

READENERGY= readenergy.o
RAND= randomnumber.o
SRC = paramod.o interface.o gendet-amlt.o matpppinit.o makeinit.o savg.o $(RAND)


#rules:
.SUFFIXES: .F90 
.SUFFIXES: .f90 
.SUFFIXES: .c

.F90.o:
	$(F90S) $(FLAGS) $(FPPFLAGS) -c -o $@ $<
.f90.o:
	$(F90S) $(FLAGS) $(FPPFLAGS) -c -o $@ $<
.f.o:
	$(F77) $(FLAGS) $(FPPFLAGS) -c -o $@ $<
.c.o:
	$(CC) -O3 -c $<

#
all: rbm readenergy

rand: $(RAND)
	$(CC) -O3 -c $(RAND)

readenergy: $(READENERGY)
	$(F90S) $(FLAGS) $(READENERGY) $(MATHLIBS) -o readenergy.x	

rbm:$(SRC)
	$(F90S) $(FLAGS) $(SRC) $(MATHLIBS) -o rbm.x

clean:
	rm -f *.o *.x *.dat fort.* *.out *.mod
