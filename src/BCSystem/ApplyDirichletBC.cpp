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
// apply dirichlet boundary condition in AsFem

#include "BCSystem/BCSystem.h"

void BCSystem::ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,Vec &U)
{

    int e,rankne,eStart,eEnd;

    int nDofsPerElmt,i,j,nDofsPerNode,iInd;
    int elDofsConn[270]={0};
    //*****************************************
    string sidename;
    double value;

    nDofsPerNode=dofHandler.GetDofsPerNode();


    for(i=0;i<bcInfo.GetBCBlockNum();i++)
    {
        if(bcInfo.GetIthBCKernelName(i+1)=="dirichlet")
        {
            MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
            MPI_Comm_size(PETSC_COMM_WORLD,&size);
            value=bcInfo.GetIthBCKernelValue(i+1);
            sidename=bcInfo.GetIthBCKernelSideName(i+1);

            iInd=bcInfo.GetIthBCKernelDofIndex(i+1);

            rankne=dofHandler.GetBCSideElmtsNum(sidename)/size;
            eStart=rank*rankne;
            eEnd=(rank+1)*rankne;

            if(rank==size-1) eEnd=dofHandler.GetBCSideElmtsNum(sidename);

            for(e=eStart;e<eEnd;e++)
            {
                dofHandler.GetLocalBCDofMap(sidename,e+1,nDofsPerElmt,elDofsConn);
                for(j=iInd;j<=nDofsPerElmt;j+=nDofsPerNode)
                {
                    VecSetValue(U,elDofsConn[j-1]-1,value,INSERT_VALUES);
                }
            }
        }
    }

    VecAssemblyBegin(U);
    VecAssemblyEnd(U);

}

