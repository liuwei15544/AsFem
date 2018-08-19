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

void BCSystem::ApplyDirichletBC(string sidename,
                                int dofindex,
                                double value,
                                Mesh &mesh,
                                DofHandler &dofHandler,
                                EquationSystem &equationSystem)
{
    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    int e,rankne,eStart,eEnd;


    int nDofsPerElmt,i,nDofsPerNode;
    int elDofsConn[270]={0};

    nDofsPerNode=dofHandler.GetDofsPerNode();
    rankne=dofHandler.GetBCSideElmtsNum(sidename)/size;
    eStart=rank*rankne;
    eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=dofHandler.GetBCSideElmtsNum(sidename);


    for(e=eStart;e<eEnd;e++)
    {
        dofHandler.GetLocalBCDofMap(sidename,e+1,nDofsPerElmt,elDofsConn);
        for(i=dofindex;i<=nDofsPerElmt;i+=nDofsPerNode)// start from 1!!!
        {
            VecSetValue(equationSystem.U,elDofsConn[i-1]-1,value,INSERT_VALUES);
        }
    }
}

