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
// Created by walkandthinker on 25.08.18.
// apply random initial value with noise in AsFem

#include "ICSystem/ICSystem.h"

void ICSystem::RandomNoiseIC(ICBlockInfo &icBlock, Mesh &mesh, Vec &U)
{
    // params=value noisevalue
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    PetscRandomCreate(PETSC_COMM_WORLD,&rnd);


    double xmin,xmax,value,noise;
    if(icBlock.ICParams.size()<2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: RandomNoiseIC need 2 parameters!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    value=icBlock.ICParams[0];
    noise=icBlock.ICParams[1];

    xmin=value-noise;
    xmax=value+noise;

    PetscRandomSetInterval(rnd,xmin,xmax);


    int i,j;
    int iStart,iEnd,rankn;
    int iInd;



    iInd=icBlock.ICBlockDofsIndex;

    rankn=mesh.GetNodesNum()/size;
    iStart=rank*rankn;
    iEnd=(rank+1)*rankn;

    if(rank==size-1) iEnd=mesh.GetNodesNum();

    for(i=iStart;i<iEnd;i++)
    {
        j=i*nDofsPerNode+iInd-1;
        PetscRandomGetValue(rnd,&value);
        VecSetValue(U,j,value,ADD_VALUES);
    }
}