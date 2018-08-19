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
// read [dofs]/[] block form input file

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadDofsName(EquationSystem &equationSystem)
{
    //*****************************
    //*** Format:
    //***     [dofs]
    //***       name=ux uy
    //***     []
    //*****************************
    string str,line,line0="dofs";
    vector<string> namelist;
    vector<int> indexlist;
    int linenum,blockstartlinenum,i;


    linenum=0;
    namelist.clear();
    indexlist.clear();
    if(IsBracketMatch(in,line0,blockstartlinenum))
    {
        GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
        getline(in,line);linenum+=1;// read [dofs]

        // remove empty lines
        getline(in,line0);linenum+=1;
        line=RemoveSpace(line0);
        while(line.size()<1)
        {
            getline(in,line0);linenum+=1;
            line=RemoveSpace(line0);
        }
        if(line.find("name=")!=string::npos)
        {
            i=line0.find("=");
            unsigned int j=i+1;
            str=line0.substr(j,line0.size()-i+1);
            namelist=SplitStr(str,' ');
            bool IsNoEmptyElement=false;
            bool FoundEmpty=false;
            while(!IsNoEmptyElement)
            {
                IsNoEmptyElement=false;
                FoundEmpty=false;
                for(j=0;j<namelist.size();j++)
                {
                    if(namelist[j].size()<1)
                    {
                        FoundEmpty=true;
                        namelist.erase(namelist.begin()+j);
                    }
                }
                if(j==namelist.size())
                {
                    if(!FoundEmpty)
                    {
                        IsNoEmptyElement=true;
                        break;
                    }
                }
            }
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'name=' can't found in line=%-3d ***\n",linenum);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name=ux uy..   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            return false;
        }
    }
    if(namelist.size()<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dofs name in line=%-3d***\n",linenum);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name=ux uy..   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        return false;
    }
    else
    {
        indexlist.clear();
        for(i=0;i<int(namelist.size());i++)
        {
            indexlist.push_back(i+1);
        }

        equationSystem.SetSolutionNameFromVector(namelist,indexlist);
        return true;
    }

}


