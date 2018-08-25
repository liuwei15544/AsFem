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
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Start initializing FE system...        ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing input system...   ***\n");
    inputSystem.InitInputSystem(args,argv);
    inputSystem.ReadInputFile(mesh,equationSystem,bcSystem,elementSystem.kernelBlockInfo,bcBlockList);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   input system initialized!            ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing dofhandler...     ***\n");
    dofHandler.CreateLocalToGlobalDofMap(mesh,equationSystem.GetDofsNumPerNode());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   dofhandler initialized!              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");



    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing equation system...***\n");
    equationSystem.SetDofsNum(dofHandler.GetDofsNum());
    equationSystem.Init();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   equation system initialized!         ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing bc system...      ***\n");
    bcSystem.InitFromBCBlockList(bcBlockList);
    bcSystem.SetDims(mesh.GetDims());
    bcSystem.SetUpBCSystem(equationSystem);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   bc system initialized!               ***\n");



    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing uel system...     ***\n");
    elementSystem.Init();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   uel system initialized!              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");

    linearSolver.InitSolver();
    nonlinearSolver.Init();


}

