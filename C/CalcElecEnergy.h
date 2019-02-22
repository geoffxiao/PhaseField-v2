#ifndef CALC_ELEC_ENERGY_H
#define CALC_ELEC_ENERGY_H

#include <pthread.h>
#include <stdio.h>
#include <string.h>

#include "Setup.h"
#include "ElectricSetup.h"
#include "HelperFunctions.h"
#include "PVector.h"
#include "ElectricWorkspace.h"
#include "TotalState.h"

void CalcElecEnergy(struct PVector * P, 
	struct ConstantsStruct * C, 
	struct ElectricSolutionStruct * S,
	struct ElectricWorkspaceStruct * E,
	struct TotalStateStruct * T);

#endif