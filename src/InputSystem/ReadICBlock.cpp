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
    // [bcs]
    //   [left]
    //     type=dirichlet
    //     dof=c1
    //     value=1.0
    //     side=left
    //   [end]
    //   [right]
    //     type=neumann
    //     dof=c2
    //     value=2.0
    //     side=right
    //   [end]
    // []

    icBlockList.clear();
    ICBlockInfo icBlockInfo;

    int linenum=0;
    string str,line,line0="ics";
    int blockstartlinenum;
    vector<double> numbers;
    bool ReadKernelBlockSuccess=false;

    bool HasDof=false,HasType=false;
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
                    bcBlockInfo.Clear();
                    bcBlockInfo.BCBlockName=str;

                    GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
                    getline(in,line0);linenum+=1;// read [kernel1]


                    HasType=false;
                    HasDof=false;
                    HasValue=false;
                    HasSideName=false;

                    // read type
                    getline(in,line0);linenum+=1;
                    line=RemoveSpace(line0);
                    while(line.compare(0,5,"[end]")!=0)
                    {
                        if(line.compare(0,5,"type=")==0)
                        {
                            HasType=true;
                            bcBlockInfo.BCBlockKernelName=line.substr(5,line.size()-5+1);
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
                            bcBlockInfo.BCBlockDofsName=line.substr(4,line.size()-4+1);
                        }
                        else if(line.compare(0,5,"side=")==0)
                        {
                            if(line.size()<8)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find side name in line-%2d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give side=side_name! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscFinalize();
                                abort();
                            }
                            HasSideName=true;
                            bcBlockInfo.sidename=line.substr(5,line.size()-5+1);
                        }
                        else if(line.compare(0,6,"value=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find value in line-%2d     ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give value=value!!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                PetscFinalize();
                                abort();
                            }
                            bcBlockInfo.BCValue=numbers[0];
                            HasValue=true;
                        }

                        getline(in,line0);linenum+=1;
                        line=RemoveSpace(line0);
                    }

                    if(!HasType)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find bc kernel type!!!    ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give bc kernel name! ***\n");
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
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find value= info !!!      ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give value=values !! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscFinalize();
                        abort();
                    }
                    else if(!HasSideName)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find side= info !!!       ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should give side=side_name! ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        PetscFinalize();
                        abort();
                    }
                    else
                    {
                        bcBlockList.push_back(bcBlockInfo);
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
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---cant find matched [kernels][] block--***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- please check your input file!!!   --***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return ReadKernelBlockSuccess;
}