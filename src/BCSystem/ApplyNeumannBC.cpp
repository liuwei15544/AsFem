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
// define neumann boundary condition in AsFem

#include "BCSystem/BCSystem.h"

#include "GaussPoint/GaussPoint.h"
#include "ShapeFuns/ShapeFuns.h"

void BCSystem::ApplyNeumannBC(Mesh &mesh,
                              DofHandler &dofHandler,
                              Vec &RHS)
{
    int e,rankne,eStart,eEnd;

    int nDofsPerElmt,i,j,nDofsPerNode,iInd;
    int elDofsConn[270]={0};
    bool HasNeumannBC=false;
    double elCoords[27][4];
    double shp[27][4],gs[125][4],xsj,JxW;
    double xi,eta,zeta;
    int ngp,Lint;
    //*****************************************
    string sidename;
    double value;

    nDofsPerNode=dofHandler.GetDofsPerNode();


    if(bcInfo.GetBCBlockNum()<1) return;

    HasNeumannBC=false;
    for(i=0;i<bcInfo.GetBCBlockNum();i++)
    {
        if(bcInfo.GetIthBCKernelName(i+1)=="neumann")
        {
            HasNeumannBC=true;
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
                mesh.GetLocalBCCoords(sidename,e+1,elCoords);
                if(nDims==1)
                {
                    // For 1D case, surface element is just node, so
                    // just directly add value to rhs
                    for(j=iInd;j<=nDofsPerElmt;j+=nDofsPerNode)
                    {
                        VecSetValue(RHS,elDofsConn[j-1]-1,value,ADD_VALUES);
                    }
                }
                else if(nDims==2)
                {
                    ngp=3;
                    // For 2D case, do 1D line integration
                    Int1D(ngp,Lint,gs);

                }
            }
        }
    }

    if(HasNeumannBC)
    {
        VecAssemblyBegin(RHS);
        VecAssemblyEnd(RHS);
    }
}
