.PHONY : All Debug Run

All : 
	gcc -std=c99 -o Main.exe -Wno-implicit-function-declaration \
	Setup.c HelperFunctions.c Main.c \
	PVector.c ElectricWorkspace.c TotalState.c \
	Eigenstrain.c ElasticWorkspace.c TimeStep.c CalcLandauEnergy.c \
	ElasticSetup.c ElectricSetup.c TimeSetup.c Energies.c \
	CalcElecEnergy.c CalcElasticEnergy.c \
	-L './lapack' -llapack -lblas -lgfortran \
	-I$(FFTW3INCLUDE) -L$(FFTW3DIR) -lfftw3 \
	-lm -lpthread
# -I '/mnt/HA/opt/fftw3/gcc/64/3.3.3/include' -L '/mnt/HA/opt/fftw3/gcc/64/3.3.3/lib' -lfftw3 -lm -lpthread \

Run :
	./Main.exe

Debug : 
	gcc -std=c99 -g -Wall -o Main.debug \
	Setup.c HelperFunctions.c Main.c \
	PVector.c ElectricWorkspace.c TotalState.c \
	Eigenstrain.c ElasticWorkspace.c TimeStep.c CalcLandauEnergy.c \
	ElasticSetup.c ElectricSetup.c TimeSetup.c Energies.c \
	CalcElecEnergy.c CalcElasticEnergy.c \
	-L './lapack' -llapack -lblas -lgfortran -lm \
	-I './fftw3/include/' -L './fftw3/lib/' -lfftw3 -lm -lpthread
	valgrind --leak-check=full --track-origins=yes ./Main.debug 1 1

Clean :
	rm -rf Main.exe Main.debug
