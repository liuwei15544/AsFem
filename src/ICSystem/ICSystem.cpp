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
// Created by walkandthinker on 25.08.18.
// define the initial condition system in AsFem

#include "ICSystem/ICSystem.h"

ICSystem::ICSystem()
{
    IsInit=false;
    icInfo.Init();
}
//*******************
void ICSystem::Init()
{
    IsInit=false;
    icInfo.Init();
}
//**********************
void ICSystem::AddICKernelBlock(ICBlockInfo &icBlockInfo)
{
    icInfo.AddICKernelBlock(icBlockInfo);
}
//*******************************
void ICSystem::InitFromICBlockList(vector<ICBlockInfo> &icBlockList)
{
    icInfo.AddICKernelBlockFromList(icBlockList);
    IsInit=true;
}
//*****************************

void ICSystem::SetUpICSystem(EquationSystem &equationSystem)
{
    icInfo.GenerateICKernelDofMap(equationSystem);
}

//*****************************
void ICSystem::PrintBCInfo() const
{
    if(IsInit)
    {
        icInfo.PrintICInfo();
    }
}
