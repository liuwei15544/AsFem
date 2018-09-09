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
// Created by walkandthinker on 09.09.18.
// compute the deformation gradient based stress material

#include "ElementSystem/ElementSystem.h"

void ElementSystem::DeformGradMaterial(const int &nDim,
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

