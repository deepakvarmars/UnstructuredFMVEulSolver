
FC = gfortran
FFLAGS = -O3  -ffree-form -Wall -Wno-conversion

TARGETS = fveul

ALL: $(TARGETS)

SRC = $(wildcard *.f95)
OBJ = $(patsubst %.f95,%.o,$(SRC))

fveul: $(OBJ)
	$(FC) -o fveul $(OBJ)

comvar.mod: comvar.f95
	$(FC) -c $(FFLAGS) comvar.f95

%.o: %.f95 comvar.mod
	$(FC) -c $(FFLAGS) $*.f95

userc.mod: comvar.f95
	$(FC) -c $(FFLAGS) comvar.f95

%.o: %.f95 userc.mod
	$(FC) -c $(FFLAGS) $*.f95

clean:
	rm -f $(OBJ) $(TARGETS) *.mod *.pdf 

clean2:
	rm -f $(OBJ) $(TARGETS) *.mod *.plt *.dat *.txt
