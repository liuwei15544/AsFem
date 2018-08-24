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
// define the FEsystem in AsFem

#include "FESystem/FESystem.h"

void FESystem::SetupFESystem()
{

}

//*********************************
void FESystem::SetJobType()
{
    if(feSystemInfo.jobtype=="static")
    {
        jobType=JobType::Static;
    }
    else if(feSystemInfo.jobtype=="transient")
    {
        jobType=JobType::Transient;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported job type!!!         ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        type=static or transient is ok! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
}
//*******************************
void FESystem::SetTimeStep()
{
    TotalTimeStep=feSystemInfo.totalstep;
}