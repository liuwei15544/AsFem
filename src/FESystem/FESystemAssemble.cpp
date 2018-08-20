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

#include "FESystem/FESystem.h"

void FESystem::AssembleLocalToGlobal(const int &iState, Mat &AMATRIX, Vec &RHS, Mat &Proj)
{
    int i,j;
    if(iState%3==0)
    {
        for(i=0;i<nDofsPerElmt;i++)
        {
            // Assemble local residual to global one
            VecSetValue(RHS,elDofsConn[i]-1,localRHS[i],ADD_VALUES);

            if(iState==6)
            {
                for(j=0;j<nDofsPerElmt;j++)
                {
                    MatSetValue(AMATRIX,elDofsConn[i]-1,elDofsConn[j]-1,localK[i][j],ADD_VALUES);
                }
            }
        }
    }
    else if(iState==8)
    {
        for(i=0;i<nNodesPerElmt;i++)
        {
            for(j=0;j<=12;j++)
            {
                MatSetValue(Proj,elConn[i],j,localProj[i][j],ADD_VALUES);
            }
        }
    }
}

//*****************************
void FESystem::FinishAssemble(const int &iState,Mat &AMATRIX, Vec &RHS, Mat &Proj)
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
        MatAssemblyBegin(Proj,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Proj,MAT_FINAL_ASSEMBLY);
    }

}


