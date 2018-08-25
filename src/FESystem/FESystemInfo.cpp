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

    jobtype="static";
    inputfilename="";
}

//***********************************
void FESystemInfo::PrintFESystemInfo() const
{
   PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: no impelementation !\n");
   PetscFinalize();
   abort();
}

