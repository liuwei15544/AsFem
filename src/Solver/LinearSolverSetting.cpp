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
// Created by walkandthinker on 21.08.18.
// settings for linear solve system

#include "Solver/LinearSolver.h"

void LinearSolver::SetKSPAbsoluteError(PetscReal atol)
{
    if(atol<1.0e-20)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:atol(ksp)=%8.2e is too small! ***\n",atol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    AbsoluteError=atol;
}

//*******************
void LinearSolver::SetKSPRelativeError(PetscReal rtol)
{
    if(rtol<1.0e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:rtol(ksp)=%8.2e is too small! ***\n",rtol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    RelativeError=rtol;
}

//****************************
void LinearSolver::SetKSPDError(PetscReal dtol)
{
    if(dtol>1.0e5||dtol<1.0e2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:dtol(ksp)=%8.2e isn't correct! ***\n",dtol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    DError=dtol;
}

//*************************
void LinearSolver::SetKSPMaxIterations(int maxiters)
{
    if(maxiters<1000||maxiters>3000000)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: maxiters=%8d isn't correct!***\n",maxiters);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    else
    {
        MaxIterations=maxiters;
    }
}



