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
