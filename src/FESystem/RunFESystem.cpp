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
// for the fe analysis procedure in AsFem

#include "FESystem/FESystem.h"

void FESystem::Run()
{
    if(feSystemInfo.jobtype=="static")
    {
        nonlinearSolver.SetPrintFlag(feSystemInfo.IsDebugOn);
        StaticAnalysis();
        if(feSystemInfo.IsProjOutput)
        {
            feSystemInfo.iState=8;
            fe.FormKR(feSystemInfo.iState,feSystemInfo.current_dt,feSystemInfo.current_time,feSystemInfo.ctan,mesh,dofHandler,elementSystem,equationSystem.U,equationSystem.V,equationSystem.AMATRIX,equationSystem.RHS,equationSystem.Proj);
            outputSystem.WriteUAndProjToVTUFile(mesh,equationSystem,equationSystem.U,equationSystem.Proj);
        }
        else
        {
            outputSystem.WriteUToVTUFile(mesh,equationSystem,equationSystem.U);
        }
    }
    else
    {
        TransientAnalysis();
    }

}

