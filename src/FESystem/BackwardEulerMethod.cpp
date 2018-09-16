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
// Created by walkandthinker on 18.08.18.
// define the backward euler time integration method in AsFem

#include "FESystem/FESystem.h"

void FESystem::BackwardEulerMethod()
{
    feSystemInfo.ctan[0]=1.0;feSystemInfo.ctan[1]=1.0/feSystemInfo.current_dt;

    feSystemInfo.iState=6;
    icSystem.ApplyInitialCondition(mesh,equationSystem.U0);

    if(feSystemInfo.IsProjOutput)
    {
        outputSystem.WriteUAndProjToVTUFile(0,mesh,equationSystem,equationSystem.U,equationSystem.Proj);
    }
    else
    {
        outputSystem.WriteUToVTUFile(0,mesh,equationSystem,equationSystem.U0);
    }

    VecCopy(equationSystem.U0,equationSystem.U);
    for(int step=1;step<=feSystemInfo.totalstep;step++)
    {
        feSystemInfo.currentstep=step;
        VecSet(equationSystem.V,0.0);
        if(!nonlinearSolver.Solve(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver,feSystemInfo))
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: transient nonlinear solve failed***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** step=%6d, iter=%6d               ***\n",step,nonlinearSolver.GetCurrentIters());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | R0|=%11.5e, | R|=%11.5e  ***\n",nonlinearSolver.GetR0Norm(),nonlinearSolver.GetRnorm());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   |dU0|=%11.5e, |dU|=%11.5e  ***\n",nonlinearSolver.GetdU0Norm(),nonlinearSolver.GetdUNorm());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | E0|=%11.5e, | E|=%11.5e  ***\n",nonlinearSolver.GetE0Norm(),nonlinearSolver.GetENorm());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        if(feSystemInfo.IsProjOutput)
        {
            feSystemInfo.iState=8;
            fe.FormKR(feSystemInfo.iState,feSystemInfo.current_dt,feSystemInfo.current_time,feSystemInfo.ctan,mesh,dofHandler,elementSystem,equationSystem.U,equationSystem.V,equationSystem.AMATRIX,equationSystem.RHS,equationSystem.Proj);
            outputSystem.WriteUAndProjToVTUFile(step,mesh,equationSystem,equationSystem.U,equationSystem.Proj);
            feSystemInfo.iState=6;
        }
        else
        {
            outputSystem.WriteUToVTUFile(step,mesh,equationSystem,equationSystem.U);
        }

        equationSystem.UpdateU0(equationSystem.U);
        feSystemInfo.current_time+=feSystemInfo.current_dt;
    }
}