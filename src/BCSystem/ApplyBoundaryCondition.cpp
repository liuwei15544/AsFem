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
// Apply boundary condition in AsFem


#include "BCSystem/BCSystem.h"

void BCSystem::ApplyBoundaryCondition(Mesh &mesh,DofHandler &dofHandler,EquationSystem &equationSystem)
{
    string sidename;
    double value;
    int ind;

    for(int i=0;i<bcInfo.GetBCBlockNum();i++)
    {
        sidename=bcInfo.GetIthBCKernelSideName(i+1);

        value=bcInfo.GetIthBCKernelValue(i+1);
        ind=bcInfo.GetIthBCKernelDofIndex(i+1);
        cout<<"dirichlet"<<endl;
        if(bcInfo.GetIthBCKernelName(i+1)=="dirichlet")
        {
            // Apply dirichlet boundary condition
            cout<<"dirichlet"<<endl;
            ApplyDirichletBC(sidename,ind,value,mesh,dofHandler,equationSystem);
        }
        else if(SingleBCBlock.BCBlockKernelName=="neumann")
        {
            // Apply neumann boundary condition

        }
    }

    VecAssemblyBegin(equationSystem.U);
    VecAssemblyEnd(equationSystem.U);
}

