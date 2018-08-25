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
// Define the output system in AsFem

#include "OutputSystem/OutputSystem.h"

OutputSystem::OutputSystem()
{
    IsInit=false;
    InputFileName="";
}

//***********************
void OutputSystem::SetInputFileName(string filename)
{
    InputFileName=filename;
    IsInit=true;
}

