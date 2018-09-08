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
// Created by walkandthinker on 19.08.18.
// Assemble the local contribution to global FE system

#include "FE/FE.h"

void FE::AssembleLocalToGlobal(const int &iState, Mat &AMATRIX, Vec &RHS, Vec &Proj)
{
    int i,j,iInd;
    if(iState%3==0)
    {

        // Assemble local residual to global one
        VecSetValues(RHS,nDofsPerElmt,elDofsConn,localRHS,ADD_VALUES);
        if(iState==6)
        {
            MatSetValues(AMATRIX,nDofsPerElmt,elDofsConn,nDofsPerElmt,elDofsConn,localK,ADD_VALUES);
        }
    }
    else if(iState==8)
    {
        for(i=0;i<nNodesPerElmt;i++)
        {
            for(j=0;j<=12;j++)
            {
                iInd=elConn[i]*13+j;
                VecSetValue(Proj,iInd,localProj[i][j],ADD_VALUES);
            }
        }
    }
}

//*****************************
void FE::FinishAssemble(const int &iState,Mat &AMATRIX, Vec &RHS, Vec &Proj)
{
    if(iState%3==0)
    {
        VecAssemblyBegin(RHS);
        VecAssemblyEnd(RHS);
        if(iState==6)
        {
            MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(AMATRIX,MAT_FINAL_ASSEMBLY);
        }
    }
    else if(iState==8)
    {
        VecAssemblyBegin(Proj);
        VecAssemblyEnd(Proj);
    }

}

//**************************************************
void FE::ModifyLocalKandRHS(DofHandler &dofHandler)
{
    // here we modify the local K and RHS
    // for dirichlet boundary condition during N-R iteration
    // if i-th dof is dirichlet bc, then [i,:]=[:,i]=0 [i,i]=1 RHS[i]=0
    int i,j,iInd,jInd;
    int i1,j1;
    for(i=1;i<=nNodesPerElmt;i++)
    {
        iInd=elConn[i-1];
        for(j=1;j<=nDofsPerNode;j++)
        {
            if(dofHandler.GetIthNodalJthDofState(iInd+1,j)==0)
            {
                // if dirichlet bc applied
                jInd=(i-1)*nDofsPerNode+j-1;
                for(i1=0;i1<nNodesPerElmt;i1++)
                {
                    localK[i1*nDofsPerElmt+jInd]=0.0;
                    localK[jInd*nDofsPerElmt+i1]=0.0;
                }
                localK[jInd*nDofsPerElmt+jInd]=1.0;
                localRHS[jInd]=0.0;
            }
        }
    }
}