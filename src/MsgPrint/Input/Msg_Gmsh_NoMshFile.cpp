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
// msg print for error of no msh file

#include "MsgPrint/MsgPrintForInput.h"

void Msg_Gmsh_NoMshFile()
{
    cout<<"*********************************************************"<<endl;
    cout<<"*** Error: can't find any .msh file !!!               ***"<<endl;
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}

