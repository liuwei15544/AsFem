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
// Created by walkandthinker on 16.08.18.
// split number from string after pos

#include "Utils/StringUtils.h"

vector<double> SplitNumAfter(string instr,int pos)
{
    string str;
    str=instr.substr(pos,instr.size()-pos+1);
    return SplitNum(str);
}



