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
// define the boundary system in AsFem

#include "BCSystem/BCSystem.h"

BCSystem::BCSystem()
{
    IsInit=false;
    bcInfo.Init();
    //BCBlockList.clear();
}

//****************************

void BCSystem::Init()
{
    IsInit=false;
    bcInfo.Init();
    //BCBlockList.clear();
}

//****************************

void BCSystem::SetDims(int dim)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d is an invalid value!!!   ***\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    nDims=dim;
    IsInit=true;
}

//*******************************
void BCSystem::AddBCKernelBlock(BCBlockInfo &bcBlockInfo)
{
    bcInfo.AddBCKernelBlock(bcBlockInfo);
}
//*******************************
void BCSystem::InitFromBCBlockList(vector<BCBlockInfo> &bcBlockList)
{
    bcInfo.AddBCKernelBlockFromList(bcBlockList);
    IsInit=true;
}
//*****************************

void BCSystem::SetUpBCSystem(EquationSystem &equationSystem)
{
    bcInfo.GenerateBCKernelDofMap(equationSystem);
}

//*****************************
void BCSystem::PrintBCInfo() const
{
    if(IsInit)
    {
        bcInfo.PrintBCInfo();
    }
}



