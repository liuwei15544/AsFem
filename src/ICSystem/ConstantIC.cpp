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
// apply constant initial condition system in AsFem

#include "ICSystem/ICSystem.h"

void ICSystem::ConstantIC(ICBlockInfo &icBlock,Mesh &mesh, Vec &U)
{
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    double value=icBlock.ICParams[0];


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
        VecSetValue(U,j,value,ADD_VALUES);
    }
}
