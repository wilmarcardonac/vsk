FC     	= gfortran
LC     	= $(FC)
EXE    	= vsk
FITSDIR = ./lib
LIBFITS = cfitsio
INCDIR	= ./Healpix_3.30/include_vsk
IDIR	= ./Healpix_3.30/include
INFGSL  = ./include/fgsl
LIBDIR	= ./Healpix_3.30/lib_vsk	
LDIR	= ./Healpix_3.30/lib
#F_FL   	= -O3 -Wall -I$(INCDIR) -I$(IDIR) -I$(INFGSL) -DGFORTRAN -march=native -mtune=native -mfpmath=sse -fno-second-underscore -fopenmp -fPIC -g
F_FL   	= -O3 -Wall -I$(INCDIR) -I$(IDIR) -I$(INFGSL) -DGFORTRAN -fno-second-underscore -fopenmp -fPIC -g
LIB_FL 	= -L$(LIBDIR) -L$(LDIR) -L$(FITSDIR) -lfgsl -lhealpix -lchealpix -lhpxgif -l$(LIBFITS)  -Wl #, -R$(FITSDIR)
#####################
OBJ   =  arrays.o fiducial.o functions.o vsk.o

def:	$(OBJ) $(OBJNR) $(OBJODE)
	$(LC) $(F_FL) $(OBJ) $(OBJNR) $(OBJODE) -o $(EXE)  $(LIB_FL)

%.o:	%.f90
	$(FC) $(F_FL) -c $<

%.o:	%.F90
	$(FC) $(F_FL) -c $<

%.o:	%.f
	$(FC) $(F_FL) -c $<

clean :
	rm -f *.o *.mod *~ fort.* *.out $(EXE)

### put dependencies here ###

vsk.o :	arrays.o functions.o fiducial.o
