FC=gfortran
SRC=$(addprefix ./, \
	2D-spline.f90 \
	_fill_b2.f90 \
	calc_potential.f90 \
	index.f90 \
	PES_3D_subroutine_ozone_L3.f90 )

# 	C6Js.f90 \
# 	C6r.f90 \
# 	Celr.f90 \
# 	cleb.f90 \
# 	faclog.f90 \

OBJ=$(addprefix ./, \
	splin.o \
	poteval.o \
	2D-spline.o \
	_fill_b2.o \
	calc_potential.o \
	index.o \
	PES_3D_subroutine_ozone_L3.o )

#   C6Js.o \
# 	C6r.o \
# 	Celr.o \
# 	cleb.o \
# 	faclog.o \

EXE=libDawesOzone.so

FFLAGS=-O2 -fpic

main :
	        rm -f $(OBJ) *.o *.mod
# 	        $(FC) -c dynamic_parameters_3D.f90 -fpic -g -ffree-line-length-0
	        $(FC) -c splin.f90 -fpic -g -ffree-line-length-0
	        $(FC) -c poteval.f90 -fpic -g -ffree-line-length-0
	        $(FC) -c $(SRC) -fpic -ffree-line-length-0
	        $(FC) -shared -fpic $(FFLAGS) $(OBJ) $(INCLUDE_PATH) $(LIB_PATH) $(LIBSHDF) $(LIBS) -o $(EXE)
	        cp $(EXE) ../$(EXE)

# $(OBJ) : Makefile
#
# %.o:%.f90
# 	        $(FC) $(FFLAGS) $(DEFS) $(INCLUDE_PATH) -c $< -o $@
# %.o:%.F90
# 	        $(FC) $(FFLAGS) $(DEFS) $(INCLUDE_PATH) -c $< -o $@

clean:
	rm -rf $(OBJ) *.o *.mod
