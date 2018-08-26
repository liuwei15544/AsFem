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
// Created by walkandthinker on 19.08.18.
// reader for initial condition block information

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadICBlock(vector<ICBlockInfo> &icBlockList)
{
    // format:
    // [ics]
    //   [ic1]
    //     type=constic
    //     dof=c1
    //     value=1.0
    //   [end]
    //   [right]
    //     type=randomic
    //     dof=c2
    //     params=valmin valmax
    //   [end]
    // []

    icBlockList.clear();
    ICBlockInfo icBlockInfo;

    int linenum=0;
    string str,line,line0="ics";
    int blockstartlinenum;
    vector<double> numbers;
    bool ReadKernelBlockSuccess=false;

    bool HasDof=false,HasType=false,HasValue=false;
    // Read the first comment line
    getline(in,line);linenum+=1;

    if(IsBracketMatch(in,line0,blockstartlinenum))
    {
        // goes inside [bcs]/[] block pair
        GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
        getline(in,line);linenum+=1;// read [bcs]


        getline(in,line);linenum+=1;// inside [kernel1]


        str=SplitStrFromBracket(line);


        while(int(line.size())>=0)
        {
            if(line.size()>0 && str.size()>0)
            {
                if(IsSubBracketMatch(in,str,blockstartlinenum))
                {
                    icBlockInfo.Clear();
                    icBlockInfo.ICBlockName=str;

                    GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
                    getline(in,line0);linenum+=1;// read [kernel1]


                    HasType=false;
                    HasDof=false;

                    // read type
                    getline(in,line0);linenum+=1;
                    line=RemoveSpace(line0);
                    while(line.compare(0,5,"[end]")!=0)
                    {
                        if(line.compare(0,5,"type=")==0)
                        {
                            HasType=true;
                            icBlockInfo.ICBlockKernelName=line.substr(5,line.size()-5+1);
                        }
                        else if(line.compare(0,4,"dof=")==0)
                        {
                            if(line.size()<5)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dof name in line-%2d   ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give dof=dof_name!!! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscFinalize();
                                abort();
                            }
                            HasDof=true;
                            icBlockInfo.ICBlockDofsName=line.substr(4,line.size()-4+1);
                        }
                        else if(line.compare(0,7,"params=")==0)
                        {
                            int i=line0.find("=")+1;
                            icBlockInfo.ICParams=SplitNum(line0.substr(i,line0.size()-i+1));
                            if(icBlockInfo.ICParams.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find value in line-%2d     ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give params=value!!! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscFinalize();
                                abort();
                            }
                            HasValue=true;
                        }

                        getline(in,line0);linenum+=1;
                        line=RemoveSpace(line0);
                    }

                    if(!HasType)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find ic kernel type!!!    ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give ic kernel name! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscFinalize();
                        abort();
                    }
                    else if(!HasDof)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dof= info !!!        ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give dof=dof_name !! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscFinalize();
                        abort();
                    }
                    else if(!HasValue)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find params= info !!!     ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give params=v1 v2..! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscFinalize();
                        abort();
                    }
                    else
                    {
                        icBlockList.push_back(icBlockInfo);
                        ReadKernelBlockSuccess=true;
                        getline(in,line);linenum+=1;// inside [materials]
                        str=RemoveSpace(line);
                        line=str;
                        if(line.size()<1)
                        {
                            while(line.size()<1)
                            {
                                getline(in,line);linenum+=1;
                                str=RemoveSpace(line);
                                line=str;
                            }
                            continue;
                        }
                        else
                        {
                            if(line.size()==2)
                            {
                                if(line.compare(0,2,"[]")==0)
                                {
                                    break;
                                }
                            }
                            continue;
                        }
                    }
                }

                getline(in,line);linenum+=1;// inside [materials]
                str=RemoveSpace(line);
                line=str;
                str=SplitStrFromBracket(line);
                if(line.compare(0,2,"[]")==0)
                {
                    break;
                }
            }
            getline(in,line);linenum+=1;// inside [materials]
            str=RemoveSpace(line);
            line=str;
            str=SplitStrFromBracket(line);
            if(line.compare(0,2,"[]")==0 && str.size()<1)
            {
                break;
            }
        }
    }

    return ReadKernelBlockSuccess;
}