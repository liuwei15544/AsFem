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
    StaticAnalysis();

    outputSystem.WriteUToVTUFile(mesh,equationSystem,equationSystem.U);
}

