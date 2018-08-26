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
// Created by walkandthinker on 24.08.18.
// define the static analysis in AsFem

#include "FESystem/FESystem.h"

void FESystem::StaticAnalysis()
{
    feSystemInfo.current_dt=1.0;
    feSystemInfo.old_dt=1.0;
    feSystemInfo.ctan[0]=1.0;
    feSystemInfo.ctan[1]=1.0;

    nonlinearSolver.Solve(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver,feSystemInfo);
}

//*********************************
void FESystem::TransientAnalysis()
{
    feSystemInfo.current_time=0.0;
    feSystemInfo.current_dt=feSystemInfo.old_dt;
    feSystemInfo.currentstep=0;
    BackwardEulerMethod();
}