#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "Energies.h"
#include "Setup.h"
#include "PVector.h"
#include "TotalState.h"
#include "TimeSetup.h"


// Find time step 1 from time step 0
void TimeStep_1stOrd(struct PVector * P_prev, 
	struct EnergiesStruct * F_prev,
	struct ConstantsStruct * C, 
	struct TimeSolutionStruct * S, 
	struct PVector * P);

#endif