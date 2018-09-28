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
// Created by walkandthinker on 28.09.18.
// apply random type initial condition system in AsFem

#include "ICSystem/ICSystem.h"

void ICSystem::RingIC(ICBlockInfo &icBlock, Mesh &mesh, Vec &U)
{
    // params=empty[0,1] or xmin xmax
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    PetscRandomCreate(PETSC_COMM_WORLD,&rnd);

    VecSet(U,0.0);

    double r1,r2,eps;
    double d1,d2,d;
    double x,y,z;
    const double sqrt2=sqrt(2.0);
    PetscScalar value;

    if(icBlock.ICParams.size()<3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 2d RingIC need 3 parameters!    ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    r1 =icBlock.ICParams[0];
    r2 =icBlock.ICParams[1];
    eps=icBlock.ICParams[2];



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

        x=mesh.GetIthNodeJthCoord(i+1,1);
        y=mesh.GetIthNodeJthCoord(i+1,2);
        z=mesh.GetIthNodeJthCoord(i+1,3);

        d=sqrt(x*x+y*y+z*z);
        d1=-(d-r1);
        d2=d-r2;

        d=max(d1,d2);

        value=-tanh(d/(sqrt2*eps));

        VecSetValue(U,j,value,INSERT_VALUES);
    }
}