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
// define the FEsystem in AsFem

#ifndef ASFEM_FESYSTEM_H
#define ASFEM_FESYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

//**********************************
//*** For AsFem's own header file
//**********************************
#include "InputSystem/InputSystem.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCSystem/BCSystem.h"
#include "ElementSystem/ElementSystem.h"
#include "FE/FE.h"
#include "EquationSystem/EquationSystem.h"

#include "FESystemInfo.h"

using namespace std;

class FESystem
{
public:
    FESystem();
    void Init(int args,char *argv[]);
    void Run();

    //****************************************
private:
    enum JobType
    {
        Static,
        Transient
    };

    enum TimeInteMethod
    {
        BackwardEuler,
        CrankNicolson,
        BDF2
    };

    enum FEAction
    {
        StaticAnalysisAction,
        TransientAnalysisAction,
        OutputAction,
        FormKRAction,
        ProjectGaussToNodalAction,
        InitHistoryAction,
        UpdateHistoryAction,
        UpdateSolutionAction
    };
    //****************************************

private:
    void SetupFESystem();
    void SetJobType();
    void SetTimeStep();
    void SetProjectionEnable();
    void SetTimeIntegrationMethod();
    void SetSolver();
    void SetNonlinearSolver();
    // TODO: add time stepping algorithm
    void SetTimeStepper();

private:
    InputSystem inputSystem;
    Mesh mesh;
    DofHandler dofHandler;
    BCSystem bcSystem;
    vector<BCBlockInfo> bcBlockList;

    ElementSystem elementSystem;
    EquationSystem equationSystem;
    FE fe;

    FESystemInfo feSystemInfo;
private:
    JobType jobType;
    TimeInteMethod timeInteMethod;
    int TotalTimeStep=0;
    double FinalTime=0.0;



private:
    bool IsFESystemInit=false;
};


#endif //ASFEM_FESYSTEM_H
