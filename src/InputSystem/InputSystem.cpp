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
// Define the input file read system in AsFem

#include "InputSystem/InputSystem.h"

InputSystem::InputSystem()
{
    IsInputSystemInit=false;
}

//*************************************

void InputSystem::InitInputSystem(int argc, char **argv)
{
    if(IsInputSystemInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: InputSystem is already init!!!  ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't init it again!!!          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    else
    {
        string str;
        if(argc>1)
        {
            str=argv[1];
            if(str.at(0)!='-'&&(str.at(1)!='i'||str.at(1)!='I'))
            {
                InputFileName=str;
            }
            else
            {
                InputFileName=str.substr(2,str.length());
            }
            in.open(InputFileName.c_str(),ios::in);
            if(!in.is_open())
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Wrong inp file name!!!                 ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            IsInputSystemInit=true;
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- Start input step -------------------***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---   Input the inp file name:");
            cin >> InputFileName;
            in.open(InputFileName.c_str(), ios::in);

            while (!in.is_open())
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- Wrong input file name!!!            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---   Input the inp file name:");
                cin >> InputFileName;
                in.open(InputFileName.c_str(), ios::in);
            }
            IsInputSystemInit=true;
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        }
    }
}
//*********************
InputSystem::InputSystem(int argc, char **argv)
{
    string str;
    if(argc>1)
    {
        str=argv[1];
        if(str.at(0)!='-'&&(str.at(1)!='i'||str.at(1)!='I'))
        {
            InputFileName=str;
        }
        else
        {
            InputFileName=str.substr(2,str.length());
        }
        in.open(InputFileName.c_str(),ios::in);
        if(!in.is_open())
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Wrong inp file name!!!                 ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
        IsInputSystemInit=true;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- Start input step -------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---   Input the inp file name:");
        cin >> InputFileName;
        in.open(InputFileName.c_str(), ios::in);

        while (!in.is_open())
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- Wrong input file name!!!            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***---   Input the inp file name:");
            cin >> InputFileName;
            in.open(InputFileName.c_str(), ios::in);
        }
        IsInputSystemInit=true;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    }
}

