#ifndef ELASTICWORKSPACE_H
#define ELASTICWORKSPACE_H

#include <complex.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct ElasticWorkspace
{
	
} ElasticWorkspace;

void FreeElasticWorkspace(struct ElasticWorkspace * U);
void InitElasticWorkspace(int Nx, int Ny, int Nz, struct ElasticWorkspace * U);

#endif