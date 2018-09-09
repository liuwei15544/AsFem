//*****************************************************
//*** AsFem: a simple finite element method program ***
//*** Copyright (C) 2018 walkandthinker             ***
//*** Contact: walkandthinker@gmail.com             ***
//*****************************************************
//* This file is part of the ASFEM framework
//* All rights reserved, see COPYRIGHT for full restrictions
//* Licensed under GPL 3.0, please see LICENSE for details
//******************************************************
//
// Created by walkandthinker on 30.08.18.
// define the linear elastic material behavior for mechanical problems in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::LinearElasticMaterial(const int &nDim,
                                          const RankTwoTensor &grad,
                                          const RankTwoTensor &strain,
                                          RankTwoTensor &stress,
                                          RankFourTensor &Jacobian)
{
    static const double E=Parameters[0];
    static const double nu=Parameters[1];
    if(Parameters.size()<2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: for mechanics, you need E,nu!!!  ***\n");
        PetscFinalize();
        abort();
    }
    Jacobian.FillFromEandNu(E,nu);
    stress=Jacobian*strain;
}