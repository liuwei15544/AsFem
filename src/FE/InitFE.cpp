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
// Initializing the FE

#include "FE/FE.h"

void FE::Init(Mesh &mesh,DofHandler &dofHandler)
{
    nDims=mesh.GetDims();
    nNodesPerElmt=mesh.GetNodesNumPerElmt();
    nDofsPerNode=dofHandler.GetDofsPerNode();
    nDofsPerElmt=nNodesPerElmt*nDofsPerNode;
    current_elmt_id=-1;
    IsInit=true;
}

//**********************************
void FE::ZeroMatAndVec(const int &iState, Mat &AMATRIX, Mat &Proj, Vec &RHS)
{
    if(iState%3==0)
    {
        VecSet(RHS,0.0);
        if(iState==6)
        {
            MatZeroEntries(AMATRIX);
        }
    }
    else if(iState==8)
    {
        MatZeroEntries(Proj);
    }
}

