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
// define the boundary system in AsFem

#include "BCSystem/BCSystem.h"

void BCSystem::ApplyConstraint(Mesh &mesh,
                               DofHandler &dofHandler,
                               Mat &AMATRIX,Vec &RHS)
{
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    int rankne,eStart,eEnd;
    int nDofsPerElmt,i,j,e,iInd,nDofsPerNode;
    int elDofsConn[270]={0};
    string sidename;
    const double PenaltyFactor=1.0e14;
    bool HasConstrainBC=false;

    if(bcInfo.GetBCBlockNum()<1) return;

    nDofsPerNode=dofHandler.GetDofsPerNode();
    HasConstrainBC=false;
    for(i=0;i<bcInfo.GetBCBlockNum();i++)
    {
        if(bcInfo.GetIthBCKernelName(i+1)=="dirichlet")
        {
            HasConstrainBC=true;
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
                    VecSetValue(RHS,elDofsConn[j-1]-1,0.0,INSERT_VALUES);
                    MatSetValue(AMATRIX,elDofsConn[j-1]-1,elDofsConn[j-1]-1,PenaltyFactor,INSERT_VALUES);
                }
            }
        }
        else if(bcInfo.GetIthBCKernelName(i+1)=="tdirichlet")
        {
            HasConstrainBC=true;
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
                    VecSetValue(RHS,elDofsConn[j-1]-1,0.0,INSERT_VALUES);
                    MatSetValue(AMATRIX,elDofsConn[j-1]-1,elDofsConn[j-1]-1,PenaltyFactor,INSERT_VALUES);
                }
            }
        }
    }

    if(HasConstrainBC)
    {
        VecAssemblyBegin(RHS);
        VecAssemblyEnd(RHS);

        MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(AMATRIX,MAT_FINAL_ASSEMBLY);
    }

}

