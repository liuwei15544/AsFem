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
// msg for no number found in mesh block

#include "MsgPrint/MsgPrintForInput.h"

void Msg_InputFile_NoNumberFound(int linenumber)
{
    printf("*** Error: can't find number in line-%-3d !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}

