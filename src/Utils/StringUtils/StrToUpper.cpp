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
// Convert str to upper case

#include <algorithm>
#include "Utils/StringUtils.h"

string StrToUpper(string instr)
{
    string outstr=instr;
    transform(outstr.begin(),outstr.end(),outstr.begin(),::toupper);
    return outstr;
}

