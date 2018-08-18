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
// split string by symbol

#include <sstream>
#include "Utils/StringUtils.h"

vector<string> SplitStr(string instr,char symbol)
{
    vector<string> outstr;
    stringstream ss(instr);
    string tok;

    outstr.clear();
    while(getline(ss,tok,symbol))
    {
        outstr.push_back(tok);
    }
    return outstr;
}
