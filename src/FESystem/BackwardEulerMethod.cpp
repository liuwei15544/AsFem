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
// define the backward euler time integration method in AsFem

#include "FESystem/FESystem.h"

void FESystem::BackwardEulerMethod()
{
    feSystemInfo.ctan[0]=1.0;feSystemInfo.ctan[1]=1.0/feSystemInfo.current_dt;

    feSystemInfo.iState=6;
    for(int step=1;step<=feSystemInfo.totalstep;step++)
    {
        VecCopy(equationSystem.U0,equationSystem.U);
        nonlinearSolver.Solve(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver,feSystemInfo);
        if(feSystemInfo.IsProjOutput)
        {
            feSystemInfo.iState=8;
            fe.FormKR(feSystemInfo.iState,feSystemInfo.current_dt,feSystemInfo.current_time,feSystemInfo.ctan,mesh,dofHandler,elementSystem,equationSystem.U,equationSystem.V,equationSystem.AMATRIX,equationSystem.RHS,equationSystem.Proj);
            outputSystem.WriteUAndProjToVTUFile(mesh,equationSystem,equationSystem.U,equationSystem.Proj);
            feSystemInfo.iState=6;
        }
        else
        {
            outputSystem.WriteUToVTUFile(mesh,equationSystem,equationSystem.U);
        }
    }
}