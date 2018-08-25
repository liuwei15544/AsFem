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
// apply circle type initial condition system in AsFem

#include "ICSystem/ICSystem.h"

void ICSystem::CircleIC(ICBlockInfo &icBlock, Mesh &mesh, Vec &U)
{
    // params=x0 y0 r v0 v1
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    if(icBlock.ICParams.size()<5)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 2d CircleIC need 5 parameters!  ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    double x0,y0,r,v0,v1;
    x0=icBlock.ICParams[0];
    y0=icBlock.ICParams[1];
    r =icBlock.ICParams[2];
    v0=icBlock.ICParams[3];
    v1=icBlock.ICParams[4];


    int i,j;
    int iStart,iEnd,rankn;
    int iInd;
    double dist,x,y;


    iInd=icBlock.ICBlockDofsIndex;

    rankn=mesh.GetNodesNum()/size;
    iStart=rank*rankn;
    iEnd=(rank+1)*rankn;

    if(rank==size-1) iEnd=mesh.GetNodesNum();

    for(i=iStart;i<iEnd;i++)
    {
        j=i*nDofsPerNode+iInd-1;
        x=mesh.GetIthNodeJthCoord(i+1,1);
        y=mesh.GetIthNodeJthCoord(i+1,2);
        dist=(x-x0)*(x-x0)+(y-y0)*(y-y0);
        dist=sqrt(dist);

        if(dist<=r)
        {
            VecSetValue(U,j,v0,ADD_VALUES);
        }
        else
        {
            VecSetValue(U,j,v0,ADD_VALUES);
        }
    }
}