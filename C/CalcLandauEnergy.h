#ifndef CALC_LANDAU_H
#define CALC_LANDAU_H

#include "TotalState.h"
#include "PVector.h"
#include "CalcLandauEnergy.h"
#include "Setup.h"

void CalcLandauEnergy(struct ConstantsStruct * C, struct TotalStateStruct * T,
	struct PVector * P, int start, int end);
void CalcLandauEnergy_MultiThreaded(struct ConstantsStruct * C, struct TotalStateStruct * T,
	struct PVector * P);
	
#endif