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
        Msg_InputSystemInitWarning();
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
                Msg_InputFileNameWrong(InputFileName);
                Msg_ExitProgram();
            }
            IsInputSystemInit=true;
        }
        else
        {
            cout<<"***---------------------------------------------------***"<<endl;
            cout<<"***--- Start the input file reading ... --------------***"<<endl;
            cout<<"***---   Enter the inp file name:";
            cin >> InputFileName;
            in.open(InputFileName.c_str(), ios::in);

            while (!in.is_open())
            {
                Msg_InputFileNameWrong(InputFileName);
                cout<<"***---   Enter the correct inp file name:";
                cin >> InputFileName;
                in.open(InputFileName.c_str(), ios::in);
            }
            IsInputSystemInit=true;
            cout<<"***---------------------------------------------------***"<<endl;

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
            Msg_InputFileNameWrong(InputFileName);
            Msg_ExitProgram();
        }
        IsInputSystemInit=true;
    }
    else
    {
        cout<<"***---------------------------------------------------***"<<endl;
        cout<<"***--- Start the input file reading ... --------------***"<<endl;
        cout<<"***---   Enter the inp file name:";
        cin >> InputFileName;
        in.open(InputFileName.c_str(), ios::in);


        while (!in.is_open())
        {
            Msg_InputFileNameWrong(InputFileName);
            cout<<"***---   Enter the correct inp file name:";
            cin >> InputFileName;
            in.open(InputFileName.c_str(), ios::in);
        }
        IsInputSystemInit=true;
        cout<<"***---------------------------------------------------***"<<endl;


    }
}
