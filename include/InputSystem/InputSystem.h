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
// Created by walkandthinker on 16.08.18.
// Define the input file read system in AsFem
//

#ifndef ASFEM_INPUTSYSTEM_H
#define ASFEM_INPUTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <fstream>

//****************************************
//*** For AsFem's own header file
//****************************************
#include "MsgPrint/MsgPrintForInput.h"
#include "MsgPrint/MsgPrintForProgram.h"

#include "Utils/StringUtils.h"


class InputSystem
{
public:
    InputSystem();
    InputSystem(int argc,char *argv[]);

    void InitInputSystem(int argc,char *argv[]);

    

private:
    bool IsInputSystemInit=false;
    string InputFileName;
    ifstream in;


};



#endif 
