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
// Created by walkandthinker on 19.04.19.
// msg print for error of no msh file

#include "MsgPrint/MsgPrintForInput.h"

void Msg_Gmsh_ImportFailed()
{
    cout<<"*********************************************************"<<endl;
    cout<<"*** Error: mesh class read gmsh file failed !!!       ***"<<endl;
    cout<<"***        please check the name of *.msh file !!!    ***"<<endl;
    cout<<"*********************************************************"<<endl;
}