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
// Convert str to low case vector

#include <algorithm>
#include "Utils/StringUtils.h"

vector<string> StrVecToLower(vector<string> instrvec)
{
    vector<string> outstrvec=instrvec;
    for(unsigned int i=0;i<outstrvec.size();i++)
    {
        string str=outstrvec[i];
        transform(str.begin(),str.end(),str.begin(),::tolower);
        outstrvec[i]=str;
    }
    return outstrvec;
}

