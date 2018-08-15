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
// Created by walkandthinker on 15.08.18.
// Generate mesh for AsFem

#include "Mesh/Mesh.h"

void Mesh::CreateMesh()
{
    if(IsMeshInfoComplete())
    {
        if(nDim==1)
        {
            Create1DMesh();
        }
        else if(nDim==2)
        {
            Create2DMesh();
        }
        else if(nDim==3)
        {
            Create3DMesh();
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't create mesh!                   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        mesh information is not complete!!!  ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }
}

