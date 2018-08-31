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
// read the system control info from input file

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadFESystemInfo(FESystemInfo &feSystemInfo)
{
    feSystemInfo.Init();

    int linenum=0;
    string str,substr,line,line0="run";
    int blockstartlinenum;
    vector<double> numbers;
    bool ReadRunBlockSuccess=false;

    bool HasType=false;
    // Read the first comment line
    getline(in,line);linenum+=1;



    if(IsBracketMatch(in,line0,blockstartlinenum))
    {
        feSystemInfo.jobtype="";
        feSystemInfo.IsDebugOn=false;
        feSystemInfo.IsProjOutput=false;

        HasType=false;
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
                    HasType=false;
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find job type !!!         ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        type=static is required !!!     ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                substr=str.substr(5,str.size()-5+1);
                if(substr.size()<5)
                {
                    HasType=false;
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported job type!!!         ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        type=static is required !!!     ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                if(substr.compare("static")==0)
                {
                    feSystemInfo.jobtype="static";
                    HasType=true;
                }
                else if(substr.compare("transient")==0)
                {
                    feSystemInfo.jobtype="transient";
                    HasType=true;
                }
            }
            else if(str.compare(0,5,"proj=")==0)
            {
                if(str.size()<9)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find proj= options        ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        proj=true[false] is required!!! ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                substr=str.substr(5,str.size()-5+1);
                if(substr.compare("true")==0)
                {
                    feSystemInfo.IsProjOutput=true;
                }
                else if(substr.compare("false")==0)
                {
                    feSystemInfo.IsProjOutput=false;
                }
                else
                {
                    feSystemInfo.IsProjOutput=false;
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported proj= options       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        proj=true[false] is required!!! ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
            }
            else if(str.compare(0,3,"dt=")==0)
            {
                if(str.size()<4)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dt= value  !!!       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dt=val1... is required !!!      ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                int i=line.find("=")+1;
                numbers=SplitNum(line.substr(i,line.size()-i+1));
                if(numbers.size()<1)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dt= value  !!!       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dt=val1... is required !!!      ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                if(numbers[0]<1.0e-13)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dt=%14.6e is two small  !!!     ***\n",numbers[0]);
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dt>=1.0e-13... is required !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                feSystemInfo.old_dt=numbers[0];
            }
            else if(str.compare(0,5,"step=")==0)
            {
                if(str.size()<6)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find step= value  !!!     ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        step=val1... is required !!!    ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                int i=line.find("=")+1;
                numbers=SplitNum(line.substr(i,line.size()-i+1));
                if(numbers.size()<1)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dt= value  !!!       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dt=val1... is required !!!      ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                if(numbers[0]<1)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: step=%6d is two small  !!!   ***\n",int(numbers[0]));
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        step>=1... is required !!!      ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                feSystemInfo.totalstep=int(numbers[0]);
            }
            else if(str.compare(0,6,"debug=")==0)
            {
                if(str.size()<10)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find debug= option  !!!   ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        debug=true or false is required!***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                int i=line.find("=")+1;
                substr=line.substr(i,line.size()-i+1);
                if(substr.size()<4)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find debug= option  !!!   ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        debug=true or false is required!***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }

                if(substr.compare(0,4,"true")==0)
                {
                    feSystemInfo.IsDebugOn=true;
                }
                else if(substr.compare(0,5,"false")==0)
                {
                    feSystemInfo.IsDebugOn=false;
                }
            }
            else if(str.compare(0,8,"maxiter=")==0)
            {
                int i=line.find("=")+1;
                numbers=SplitNum(line.substr(i,line.size()-i+1));
                if(numbers.size()<1)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find dt= value  !!!       ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dt=val1... is required !!!      ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                if(numbers[0]<1)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: step=%6d is two small  !!!   ***\n",int(numbers[0]));
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        step>=1... is required !!!      ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                feSystemInfo.MaxNonlinearIter=int(numbers[0]);
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

        if(HasType)
        {
            ReadRunBlockSuccess=true;
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---[run] block info is incomplete!!! ---***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- read run block failed !!!        ---***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---cant find matched [run][] block   ---***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- please check your input file!!!  ---***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return ReadRunBlockSuccess;
}

