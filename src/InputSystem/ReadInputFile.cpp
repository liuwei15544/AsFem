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
// Created by walkandthinker on 16.08.18.
// Define the input file read system in AsFem

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadInputFile(Mesh &mesh,
                                EquationSystem &equationSystem,
                                FESystemInfo &feSystemInfo,
                                KernelBlockInfo &kernelBlockInfo,
                                vector<BCBlockInfo> &bcBlockList,
                                vector<ICBlockInfo> &icBlockList)
{
    if(!ReadMeshBlock(mesh))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: read mesh block failed!!!       ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(!ReadDofsName(equationSystem))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: read dofs block failed!!!       ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(!ReadKernelBlock(kernelBlockInfo))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: read kernel block failed!!!     ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

//    if(!ReadBoundaryBlock(bcBlockList))
//    {
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: read boundary block failed!!!   ***\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
//        PetscFinalize();
//        abort();
//    }

    ReadBoundaryBlock(bcBlockList);

//    if(!ReadICBlock(icBlockList))
//    {
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: read IC block failed!!!         ***\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
//        PetscFinalize();
//        abort();
//    }

    ReadICBlock(icBlockList);

    if(!ReadFESystemInfo(feSystemInfo))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: read run block failed!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

}

