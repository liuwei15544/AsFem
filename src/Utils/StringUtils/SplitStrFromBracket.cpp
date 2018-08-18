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

string SplitStrFromBracket(string &inputstr)
{
    string str,str1;
    str1=RemoveSpace(inputstr);

    if(str1.size()<2)
    {
        str="";
    }
    else
    {
        if((str1.at(0)=='[') && (str1.at(str1.size()-1)==']') && str1.size()>2)
        {
            str=str1.substr(1,str1.size()-2);
        }
        else
        {
            str="";
        }
    }

    return str;
}

