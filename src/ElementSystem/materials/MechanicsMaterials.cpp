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
// define the material system for mechanical problems in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::MechanicsMaterials(const int &nDim,const RankTwoTensor &grad,
                                      RankTwoTensor &strain,
                                      RankTwoTensor &stress,
                                      RankFourTensor &Jacobian)
{
    ComputeStrain(nDim,grad,strain);
    switch (ActiveUmatIndex)
    {
        case linearelastic:
            LinearElasticMaterial(nDim,grad,strain,stress,Jacobian);
            break;
        case deformgrad:
            DeformGradMaterial(nDim,grad,strain,stress,Jacobian);
            break;
        case neohookean:
            NeoHookeanMaterial(nDim,grad,strain,stress,Jacobian);
            break;
        default:
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported umat for mechanics!!***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
    }
}

//************************************
void ElementSystem::ComputeStrain(const int &nDim,
                                  const RankTwoTensor &grad,
                                  RankTwoTensor &strain)
{
    RankTwoTensor F(0.0),Ft(0.0),I(0.0);

    if(strainMode==small)
    {
        F=grad.transpose();

        strain=0.5*(F+grad);
    }
    else
    {
        I.IdentityEntities();

        F=grad+I;
        Ft=F.transpose();

        strain=0.5*(Ft*F-I);
    }
}