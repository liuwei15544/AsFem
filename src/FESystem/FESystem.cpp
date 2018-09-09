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
    inputSystem.ReadInputFile(mesh,
                              equationSystem,
                              feSystemInfo,
                              elementSystem.kernelBlockInfo,
                              bcBlockList,icBlockList);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   input system initialized!            ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");



    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing dofhandler...     ***\n");
    dofHandler.CreateLocalToGlobalDofMap(mesh,equationSystem.GetDofsNumPerNode());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   dofhandler initialized!              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");


    feSystemInfo.nDofs=dofHandler.GetDofsNum();
    feSystemInfo.nNodes=mesh.GetNodesNum();
    feSystemInfo.nElmts=mesh.GetElmtsNum();

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing equation system...***\n");
    equationSystem.SetDofsNum(dofHandler.GetDofsNum());
    equationSystem.Init();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   equation system initialized!         ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");



    if(bcBlockList.size()>0)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing bc system...      ***\n");
        bcSystem.InitFromBCBlockList(bcBlockList);
        bcSystem.SetDims(mesh.GetDims());
        bcSystem.SetUpBCSystem(equationSystem);
        dofHandler.SetNodalDofActiveState(mesh,bcSystem.bcInfo);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   bc system initialized!               ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    }


    if(icBlockList.size()>0)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing ic system...      ***\n");
        icSystem.InitFromICBlockList(icBlockList);
        icSystem.SetDofsNumPerNode(dofHandler.GetDofsPerNode());
        icSystem.SetUpICSystem(equationSystem);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   ic system initialized!               ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    }



    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing uel system...     ***\n");
    elementSystem.Init();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   uel system initialized!              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing linear solver...  ***\n");
    linearSolver.InitSolver();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   linear solver initialized!           ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");



    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing nonlinear solver..***\n");
    nonlinearSolver.Init();
    nonlinearSolver.SetMaxIters(feSystemInfo.MaxNonlinearIter);
    nonlinearSolver.SetSolverType(feSystemInfo.NonLinearSolver);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nonlinear solver initialized!        ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   start initializing output system...  ***\n");
    outputSystem.SetInputFileName(inputSystem.GetInputFileName());
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   output system initialized!           ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** FE system initialized!                 ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

    feSystemInfo.PrintFESystemInfo();

}

