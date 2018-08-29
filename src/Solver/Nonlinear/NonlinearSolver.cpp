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
// define nonlinear solve for Ax=F

#include "Solver/NonlinearSolver.h"

NonlinearSolver::NonlinearSolver()
{
    SolverType=newtonraphson;
    Rtol_R = 1.0e-7; Atol_R=5.0e-10;

    Rtol_dU= 1.0e-8;Atol_dU=1.0e-11;

    Rtol_E =1.0e-12; Atol_E=1.0e-16;

    MaxIters=50;iters=0;
    IsInit=true;
    PrintIterationInfo=false;
}

void NonlinearSolver::Init()
{
    SolverType=newtonraphson;
    Rtol_R =5.0e-8; Atol_R=1.0e-9;

    Rtol_dU=1.0e-8;Atol_dU=1.0e-11;

    Rtol_E =1.0e-16;Atol_E=1.0e-18;

    MaxIters=50;iters=0;
    IsInit=true;
}


//**********************************

void NonlinearSolver::PrintNonlinearSolverInfo() const
{
    if(IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Nonlinear solver info:                 ***\n");
        if(SolverType==newtonraphson)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   method= newton-raphson(default)      ***\n");
        }
        else if(SolverType==modifynewtonraphson)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   method= modified-newton              ***\n");
        }
        else if(SolverType==arclength)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   method= arc-length method            ***\n");
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   max iteration=%3d                    ***\n",GetMaxIters());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   relative error of  R=%14.5e  ***\n",GetRtolOfResidual());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   absolute error of  R=%14.5e  ***\n",GetAtolOfResidual());

        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   relative error of dU=%14.5e  ***\n",GetRtolOfdU());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   absolute error of dU=%14.5e  ***\n",GetAtolOfdU());

        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   relative error of  E=%14.5e  ***\n",GetRtolOfEnergy());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   absolute error of  E=%14.5e  ***\n",GetAtolOfEnergy());

    }
}


