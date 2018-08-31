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
// Created by walkandthinker on 24.08.18.
// define the static analysis in AsFem

#include "FESystem/FESystem.h"

void FESystem::StaticAnalysis()
{
    feSystemInfo.current_dt=1.0;
    feSystemInfo.old_dt=1.0;
    feSystemInfo.ctan[0]=1.0;
    feSystemInfo.ctan[1]=1.0;

    if(!nonlinearSolver.Solve(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver,feSystemInfo))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: static nonlinear solve failed!!!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** iter=%6d                            ***\n",nonlinearSolver.GetCurrentIters());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | R0|=%11.5e, | R|=%11.5e  ***\n",nonlinearSolver.GetR0Norm(),nonlinearSolver.GetRnorm());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   |dU0|=%11.5e, |dU|=%11.5e  ***\n",nonlinearSolver.GetdU0Norm(),nonlinearSolver.GetdUNorm());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | E0|=%11.5e, | E|=%11.5e  ***\n",nonlinearSolver.GetE0Norm(),nonlinearSolver.GetENorm());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
}

//*********************************
void FESystem::TransientAnalysis()
{
    feSystemInfo.current_time=0.0;
    feSystemInfo.current_dt=feSystemInfo.old_dt;
    feSystemInfo.currentstep=0;
    BackwardEulerMethod();
}