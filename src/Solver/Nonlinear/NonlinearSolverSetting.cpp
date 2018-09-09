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
// settings for nonlinear solver

#include "Solver/NonlinearSolver.h"

//********************************
void NonlinearSolver::SetMaxIters(int maxiters)
{
    if(maxiters<2 || maxiters>10000)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: maxiters=%5d is invalid!!!    ***\n",maxiters);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    MaxIters=maxiters;
}

//*********************************************
void NonlinearSolver::SetSolverType(string solvertype)
{
    if(solvertype=="newtonraphson")
    {
        SolverType=newtonraphson;
    }
    else if(solvertype=="modifynewtonraphson")
    {
        SolverType=modifynewtonraphson;
    }
    else if(solvertype=="arclength")
    {
        SolverType=arclength;
    }
    else if(solvertype=="linesearch")
    {
        SolverType=netownwithlinesearch;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unknown nonlinear solver type!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
}

//**************************
void NonlinearSolver::SetRtolOfResidual(double rtol)
{
    if(rtol<1.0e-12)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:rtol(R)=%10.2e is too small! ***\n",rtol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    Rtol_R=rtol;
}
//***************
void NonlinearSolver::SetAtolOfResidual(double atol)
{
    if(atol<1.0e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:atol(R)=%10.2e is too small! ***\n",atol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    Atol_R=atol;
}

//**********************
void NonlinearSolver::SetRtolOfdU(double rtol)
{
    if(rtol<1.0e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:rtol(dU)=%10.2e is too small!***\n",rtol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    Rtol_dU=rtol;
}
//*********************
void NonlinearSolver::SetAtolOfdU(double atol)
{
    if(atol<1.0e-15)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:atol(dU)=%10.2e is too small!***\n",atol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    Atol_dU=atol;
}
//********************************
void NonlinearSolver::SetRtolOfEnergy(double rtol)
{
    if(rtol<1.0e-12)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:rtol(E)=%10.2e is too small! ***\n",rtol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    Rtol_E=rtol;
}

//**********************
void NonlinearSolver::SetAtolOfEnergy(double atol)
{
    if(atol<1.0e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:atol(E)=%10.2e is too small! ***\n",atol);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    Atol_E=atol;
}



