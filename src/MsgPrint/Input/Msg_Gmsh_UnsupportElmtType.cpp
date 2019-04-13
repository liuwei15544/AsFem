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
// Define the msg printer for unsupported gmsh elmt type

#include "MsgPrint/MsgPrintForInput.h"

void Msg_Gmsh_UnsupportElmtType(int type)
{
    cout<<"*********************************************************"<<endl;
    cout<<"*** Error: unsupported gmsh elmt type(="<<setw(6)<<type;
    cout<<" )       ***"<<endl;
    cout<<"*********************************************************"<<endl;
}




