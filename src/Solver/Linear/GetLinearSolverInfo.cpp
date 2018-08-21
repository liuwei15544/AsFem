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
// get info for linear solve system

#include "Solver/LinearSolver.h"

int LinearSolver::GetKSPIterations() const
{
    PetscInt CurrentIters;
    if(CurrentSolveFinished)
    {
        KSPGetIterationNumber(ksp,&CurrentIters);
        return CurrentIters;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: current ksp solve isn't finished***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get iteration information ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return MaxIterations;
}

//****************
int LinearSolver::GetKSPMaxIterations() const
{
    return MaxIterations;
}

//*****************
double LinearSolver::GetKSPAbsoluteError() const
{
    return AbsoluteError;
}
//*****************
double LinearSolver::GetKSPRelativeError() const
{
    return RelativeError;
}
//*****************
double LinearSolver::GetKSPDivergenceError() const
{
    return DError;
}