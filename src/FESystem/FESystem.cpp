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
// define the FEsystem in AsFem

#include "FESystem/FESystem.h"

FESystem::FESystem()
{
    IsFESystemInit=false;
}

//*********************************
void FESystem::Init(int args, char **argv)
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Start initializing input system...     ***\n");

    inputSystem.InitInputSystem(args,argv);
    inputSystem.ReadInputFile(mesh,equationSystem,bcSystem,elementSystem.kernelBlockInfo,bcBlockList);

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Input system initialized!              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");

    //dofHandler.CreateLocalToGlobalDofMap(mesh,equationSystem.GetDofsNumPerNode());


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Start initializing equation system...  ***\n");
    equationSystem.SetDofsNum(dofHandler.GetDofsNum());
    equationSystem.Init();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Equation system initialized!           ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Start initializing bc system...        ***\n");
    bcSystem.InitFromBCBlockList(bcBlockList);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** BC system initialized!                 ***\n");



    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Start initializing UEL system...       ***\n");
    elementSystem.Init();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** UEL system initialized!                ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");


}

