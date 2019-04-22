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
// Define the msg printer for input file read system in AsFem

#ifndef ASFEM_MSGPRINTFORINPUT_H
#define ASFEM_MSGPRINTFORINPUT_H


#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

void Msg_InputFileNameWrong(string filename);
void Msg_InputSystemInitWarning();


// For gmsh
void Msg_Gmsh_UnsupportElmtType(int type);
void Msg_Gmsh_NoMshFile();
void Msg_Gmsh_ImportFailed();

// For mesh block msg
void Msg_InputFile_MeshBlock_NoDim(int linenumber);
void Msg_InputFile_MeshBlock_DimInvalid(int linenumber);
void Msg_InputFile_MeshBlock_XminNotFound(int linenumber);
void Msg_InputFile_MeshBlock_XmaxNotFound(int linenumber);
void Msg_InputFile_MeshBlock_YminNotFound(int linenumber);
void Msg_InputFile_MeshBlock_YmaxNotFound(int linenumber);
void Msg_InputFile_MeshBlock_ZminNotFound(int linenumber);
void Msg_InputFile_MeshBlock_ZmaxNotFound(int linenumber);

void Msg_InputFile_NoNumberFound(int linenumber);


#endif 
