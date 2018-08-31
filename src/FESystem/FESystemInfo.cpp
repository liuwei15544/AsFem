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
// define the FESystem info in AsFem

#include "FESystem/FESystemInfo.h"

FESystemInfo::FESystemInfo()
{
    current_dt=0.0;old_dt=0.0;
    current_time=0.0;
    totalstep=0;
    iState=6;

    IsProjOutput=false;

    jobtype="static";
    inputfilename="";
}

//*********************************
void FESystemInfo::Init()
{
    current_dt=0.0;old_dt=0.0;
    current_time=0.0;
    totalstep=0;
    iState=6;

    ctan[0]=1.0;ctan[1]=1.0;

    IsProjOutput=false;

    nDofs=0;nNodes=0;nElmts=0;
    MaxNonlinearIter=50;

    jobtype="static";
    inputfilename="";
}

//***********************************
void FESystemInfo::PrintFESystemInfo() const
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** FE system information:                 ***\n");
    if(jobtype=="static")
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** job type= static analysis              ***\n");
    }
    else if(jobtype=="transient")
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** job type= transient analysis           ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** total step= %6d                     ***\n",totalstep);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** delta t   =%14.6e              ***\n",old_dt);
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** number of dofs =%8d               ***\n",nDofs);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** number of nodes=%8d               ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** number of elmts=%8d               ***\n",nElmts);

    if(IsProjOutput==true)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** projection output=true                 ***\n");
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** projection output=false                ***\n");
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** maximum iters=%3d                      ***\n",MaxNonlinearIter);

    if(IsDebugOn)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** debug mode=true                        ***\n");
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** debug mode=false                       ***\n");
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}

