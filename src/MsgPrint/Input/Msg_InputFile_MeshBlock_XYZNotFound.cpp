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
// msg for no x, y or z found in mesh block

#include "MsgPrint/MsgPrintForInput.h"

void Msg_InputFile_MeshBlock_XminNotFound(int linenumber)
{
    printf("*** Error: can't find xmin=  in line-%3d  !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}

void Msg_InputFile_MeshBlock_XmaxNotFound(int linenumber)
{
    printf("*** Error: can't find xmax=  in line-%3d  !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}
void Msg_InputFile_MeshBlock_YminNotFound(int linenumber)
{
    printf("*** Error: can't find ymin=  in line-%3d  !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}
void Msg_InputFile_MeshBlock_YmaxNotFound(int linenumber)
{
    printf("*** Error: can't find ymax=  in line-%3d  !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}
void Msg_InputFile_MeshBlock_ZminNotFound(int linenumber)
{
    printf("*** Error: can't find zmin=  in line-%3d  !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}
void Msg_InputFile_MeshBlock_ZmaxNotFound(int linenumber)
{
    printf("*** Error: can't find zmax=  in line-%3d  !!!         ***\n",linenumber);
    cout<<"***        please check your input file !!!           ***"<<endl;
    cout<<"*********************************************************"<<endl;
}
