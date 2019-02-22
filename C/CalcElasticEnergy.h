#ifndef CALC_ELASTIC_H
#define CALC_ELASTIC_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Setup.h"
#include "ElasticSetup.h"
#include "HelperFunctions.h"
#include "PVector.h"
#include "TotalState.h"
#include "ElasticWorkspace.h"
#include "Eigenstrain.h"

void CalcElasticEnergy(struct ElasticWorkspace * U, 
	struct PVector * P,
	struct EigenstrainStruct * ES,
	struct ConstantsStruct * C, 
	struct ElasticSolutionStruct * S,
	struct TotalStateStruct * T);


#endif
