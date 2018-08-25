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

#ifndef ASFEM_INPUTSYSTEM_H
#define ASFEM_INPUTSYSTEM_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "petsc.h"

// For AsFem's own header file
#include "Utils/StringUtils.h"
#include "Mesh/Mesh.h"
#include "EquationSystem/EquationSystem.h"
#include "ElementSystem/KernelBlockInfo.h"
#include "BCSystem/BCSystem.h"

using namespace std;

class InputSystem
{
public:
    InputSystem();
    InputSystem(int argc,char *argv[]);

    void InitInputSystem(int argc,char *argv[]);

    bool ReadInputFile(Mesh &mesh,
                       EquationSystem &equationSystem,
                       BCSystem &bcSystem,
                       KernelBlockInfo &kernelBlockInfo,
                       vector<BCBlockInfo> &bcBlockList);

    string GetInputFileName() const { return InputFileName;}

private:
    bool ReadMeshBlock(Mesh &mesh);
    bool ReadDofsName(EquationSystem &equationSystem);
    bool ReadKernelBlock(KernelBlockInfo &kernelBlockInfo);
    bool ReadBoundaryBlock(vector<BCBlockInfo> &bcBlockList);


private:
    bool IsInputSystemInit=false;
    string InputFileName;
    ifstream in;

};


#endif //ASFEM_INPUTSYSTEM_H
