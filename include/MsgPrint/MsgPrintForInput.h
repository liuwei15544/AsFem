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

#endif 
