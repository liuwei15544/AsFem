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
// define the uel system in AsFem

#include "ElementSystem/ElementSystem.h"

ElementSystem::ElementSystem()
{
    kernelBlockInfo.Init();
    IsInit=false;
    ActiveUelIndex =-10000;
    ActiveUmatIndex=-10000;
}

//******************************
void ElementSystem::Init()
{
    SetUmatIndex();
    SetUelIndex();
    Parameters=kernelBlockInfo.params;
    IsInit=true;
}

