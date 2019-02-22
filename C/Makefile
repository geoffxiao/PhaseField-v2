.PHONY : All Debug Profile MemCheck Clean

All : 
	gcc -std=c99 -o Main.exe -Wno-implicit-function-declaration \
	Setup.c HelperFunctions.c Main.c \
	PVector.c ElectricWorkspace.c TotalState.c \
	Eigenstrain.c ElasticWorkspace.c TimeStep.c CalcLandauEnergy.c \
	ElasticSetup.c ElectricSetup.c TimeSetup.c Energies.c \
	CalcElecEnergy.c CalcElasticEnergy.c \
	-L './lapack' -llapack -lblas -lgfortran \
	-I './fftw3/include/' -L './fftw3/lib/' -lfftw3 -lm -lpthread
	./Main.exe

Debug : 
	gcc -std=c99 -pg -Wall -o Main.debug \
	Setup.c HelperFunctions.c Main.c \
	PVector.c ElectricWorkspace.c TotalState.c \
	Eigenstrain.c ElasticWorkspace.c TimeStep.c CalcLandauEnergy.c \
	ElasticSetup.c ElectricSetup.c TimeSetup.c Energies.c \
	CalcElecEnergy.c CalcElasticEnergy.c \
	-L './lapack' -llapack -lblas -lgfortran -lm \
	-I './fftw3/include/' -L './fftw3/lib/' -lfftw3 -lm -lpthread

Profile :
	./Main.debug
	gprof Main.debug gmon.out > output
	
MemCheck:
	valgrind --leak-check=full --track-origins=yes ./Main.debug

Clean :
	rm -rf Main.exe Main.debug output
