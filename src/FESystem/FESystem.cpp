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

#include "FESystem/FESystem.h"

FESystem::FESystem()
{
    IsFESystemInit=false;
}

//*********************************
void FESystem::Init(int args, char **argv)
{
    inputSystem.InitInputSystem(args,argv);
    inputSystem.ReadInputFile(mesh,equationSystem,bcSystem,elementSystem.kernelBlockInfo,bcBlockList);

    equationSystem.Init();
    bcSystem.InitFromBCBlockList(bcBlockList);
    elementSystem.Init();

    bcSystem.PrintBCInfo();
    equationSystem.PrintSolutionNameMap();
}

