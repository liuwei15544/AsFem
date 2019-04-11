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
// several utils for string operation

#include "Utils/StringUtils.h"

string RemoveSpace(string instr)
{
    if(instr.size()<=1)
    {
        return instr;
    }
    unsigned int i,length;
    char ch[1000];
    length=0;
    for(i=0;i<instr.size();i++)
    {
        if(instr.at(i)!=' ')
        {
            ch[length]=instr.at(i);
            length+=1;
        }
    }

    string outstr;
    outstr.clear();
    for(i=0;i<length;i++)
    {
        if(ch[i]=='\n')
        {
            break;
        }
        outstr+=ch[i];
    }
    return outstr;
}


