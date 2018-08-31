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
// define the free energy materials for CahnHilliard model in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::FreeEnergyMaterials(const double &conc, double &f, double &dfdc, double &d2fdc2)
{
    switch (ActiveUmatIndex)
    {
        case freeenergy:
            FreeEnergyForCahnHilliard(conc,f,dfdc,d2fdc2);
            break;
        default:
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported umat for free energy***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
    }
}