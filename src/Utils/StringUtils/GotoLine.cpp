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
//*** Split string from a bracked                           ***
//*** i.e. [test]-->result= test                            ***

#include "Utils/StringUtils.h"

void GotoLine(ifstream &in,int linenum)
{
    string line;
    in.clear();
    in.seekg(ios::beg);
    for(int i=0;i<linenum-1;i++)
    {
        getline(in,line);
    }
}

