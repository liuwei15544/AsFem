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
// Created by walkandthinker on 22.04.19.
// msg for no dim found in mesh block

#include "MsgPrint/MsgPrintForInput.h"

void Msg_InputFile_MeshBlock_NoDim(int linenumber)
{
    printf("*** Error: can't find 'dim=' line-%-3d                ***\n",linenumber);
    cout<<"***        you should given dim=1,2,3 !!!             ***"<<endl;
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}

void Msg_InputFile_MeshBlock_DimInvalid(int linenumber)
{
    printf("*** Error: dim=%2d is invalid in line-%-3d!!!         ***\n",linenumber);
    cout<<"***        you should given dim=1,2,3 !!!             ***"<<endl;
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}

