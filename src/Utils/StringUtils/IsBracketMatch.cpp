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
// Check whether the bracket pair is match or not

#include "Utils/StringUtils.h"

bool IsBracketMatch(ifstream &in,string &bracketstr,int &startline)
{
    string line,str;
    bool FoundHead=false,FoundEnd=false;
    string line0='['+bracketstr+']';
    int endline=0;
    startline=0;
    int linenum=0;
    in.clear();
    in.seekg(0, ios::beg);
    while(!in.eof())
    {
        getline(in,line);linenum+=1;
        str=RemoveSpace(line);
        if(str.compare(line0)==0)
        {
            FoundHead=true;


            startline=linenum;
            while(!in.eof())
            {
                getline(in,line);linenum+=1;
                str=RemoveSpace(line);

                if(str.compare("[]")==0)
                {
                    FoundEnd= true;
                    endline=linenum;
                    break;
                }
            }
        }
    }

    if((FoundHead && FoundEnd)&&(endline>startline))
    {
        //cout<<"start="<<startline<<":end="<<endline<<endl;
        return true;
    }
    else
    {
        return false;
    }
}