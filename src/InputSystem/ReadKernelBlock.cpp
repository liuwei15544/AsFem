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
// Created by walkandthinker on 18.08.18.
// implement kernel block reader of AsFem

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadKernelBlock(KernelBlockInfo &kernelBlockInfo)
{
    // format : [kernel]
    //            type=mechanics
    //            mate=neohookean
    //            params=1.0e5 0.3
    //          []
    kernelBlockInfo.Init();

    int linenum=0;
    string str,line,line0="kernel",substr;
    int blockstartlinenum;
    vector<double> numbers;
    bool ReadKernelBlockSuccess=false;

    bool HasKernelType=false;
    bool HasStrain=false,HasMate=false;
    // Read the first comment line
    getline(in,line);linenum+=1;



    if(IsBracketMatch(in,line0,blockstartlinenum))
    {
        kernelBlockInfo.MaterialKernelName="";
        HasKernelType=false;
        HasStrain=false;
        HasMate=false;
        // goes inside [kernels]/[] block pair
        GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
        getline(in,line);linenum+=1;// read [kernel]

        str=RemoveSpace(line);
        while (str.size()<1)
        {
            getline(in,line);linenum+=1;
            str=RemoveSpace(line);
        }
        while(str.compare(0,2,"[]")!=0)
        {
            if(str.compare(0,5,"type=")==0)
            {
                if(str.size()<8)
                {
                    HasKernelType=false;
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find element name !!!     ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        type=element name is required!!!***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                kernelBlockInfo.ElementName=str.substr(5,str.size()-5+1);
                HasKernelType=true;
            }
            else if(str.compare(0,5,"mate=")==0)
            {
                if(str.size()<8)
                {
                    HasKernelType=false;
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find mate name !!!        ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        mate=mate name is required !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                kernelBlockInfo.MaterialKernelName=str.substr(5,str.size()-5+1);
                HasMate=true;
            }
            else if(str.compare(0,7,"params=")==0)
            {
                if(str.size()<8)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find parameters !!!       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        params=val1... is required !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                int i=line.find("=")+1;
                kernelBlockInfo.params=SplitNum(line.substr(i,line.size()-i+1));
            }
            else if(str.compare(0,7,"strain=")==0)
            {
                if(str.size()<8)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find parameters !!!       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        params=val1... is required !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                int i=str.find("=")+1;
                substr=str.substr(i,str.size()-i+1);

                if(substr.compare(0,5,"small")==0)
                {
                    kernelBlockInfo.strain="small";
                    HasStrain=true;
                }
                else if(substr.compare(0,6,"finite")==0)
                {
                    kernelBlockInfo.strain="finite";
                    HasStrain=true;
                }
            }


            getline(in,line);linenum+=1;// inside [kernels]
            str=RemoveSpace(line);
            while (str.size()<1)
            {
                getline(in,line);linenum+=1;
                str=RemoveSpace(line);
            }

            if(str.compare(0,2,"[]")==0)
            {
                break;
            }
        }

        if(HasKernelType)
        {
            ReadKernelBlockSuccess=true;
            if(kernelBlockInfo.ElementName.compare("mechanics")==0)
            {
                if(!HasStrain)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: for the mechanics problem       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        strain=small/finite is required!***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    PetscFinalize();
                    abort();
                }

                if(!HasMate)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: for the mechanics problem       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        mate=linear/other is required!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    PetscFinalize();
                    abort();
                }
            }
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---kernel block info is incomplete!!! --***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- read kernel block failed !!!      --***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---cant find matched [kernel][] block --***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- please check your input file!!!   --***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return ReadKernelBlockSuccess;
}

