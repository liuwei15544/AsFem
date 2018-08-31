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
// define linear solve for Ax=F

#include "Solver/LinearSolver.h"

LinearSolver::LinearSolver()
{
    // Define the default options for ksp solver
    PreconditionType="lu";
    MaxIterations=2000000;
    AbsoluteError=1.0e-12;
    RelativeError=1.0e-7;
    MaxGMRSRestat=1201;
    DError=1.0e5;

    IsSolverInit=false;
    CurrentSolveFinished=false;
}

//*******************************+
bool LinearSolver::InitSolver()
{
    if(IsSolverInit)
    {
        //ierr=KSPSetTolerances(ksp,RelativeError,AbsoluteError,DError,MaxIterations);CHKERRQ(ierr);
        ierr=KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
        ierr=KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr=KSPGMRESSetRestart(ksp,MaxGMRSRestat);
        ierr=KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        ierr=PCSetType(pc,PCLU);CHKERRQ(ierr);
        return ierr;
    }
    else
    {
        ierr=KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
        ierr=KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
        //ierr=KSPSetTolerances(ksp,RelativeError,AbsoluteError,DError,MaxIterations);CHKERRQ(ierr);
        ierr=KSPSetFromOptions(ksp);CHKERRQ(ierr);

        ierr=KSPGMRESSetRestart(ksp,MaxGMRSRestat);
        ierr=KSPGetPC(ksp,&pc);CHKERRQ(ierr);
        ierr=PCSetType(pc,PCLU);CHKERRQ(ierr);

        IsSolverInit=true;

        return IsSolverInit;
    }
}

//*********************************************
//*** Solve Ax=F equations                  ***
//*********************************************
bool LinearSolver::Solve(Mat &A,Vec &x,Vec &F)
{
    CurrentSolveFinished=false;
    ierr=KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
    ierr=KSPSolve(ksp,F,x);CHKERRQ(ierr);
    CurrentSolveFinished=true;

    return CurrentSolveFinished;
}

//**************************************
//*** Release Memory
//**************************************
void LinearSolver::Release()
{
    if(IsSolverInit)
    {
        ierr=KSPDestroy(&ksp);
    }
}

//**************************************
//*** Print out solver information
//**************************************
void LinearSolver::PrintSolverInfo() const
{
    if(!IsSolverInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print solver information! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        ksp solver isn't initialized!!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** KSP solver information:                ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   max iterations=%10d            ***\n",GetKSPMaxIterations());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   relative tolerance  =%15.5e ***\n",GetKSPRelativeError());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   absolute tolerance  =%15.5e ***\n",GetKSPAbsoluteError());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   divergence tolerance=%15.5e ***\n",GetKSPDivergenceError());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        // TODO: add precondition information
    }
}




