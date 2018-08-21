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
// define the FEsystem in AsFem

#ifndef ASFEM_FESYSTEM_H
#define ASFEM_FESYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

//**********************************
//*** For AsFem's own header file
//**********************************
#include "InputSystem/InputSystem.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "BCSystem/BCSystem.h"
#include "FE/FE.h"

using namespace std;

class FESystem
{
public:
    FESystem();
    void Init(int args,char *argv[]);

private:
    InputSystem inputSystem;
    Mesh mesh;
    DofHandler dofHandler;
    BCSystem bcSystem;
    FE fe;
};


#endif //ASFEM_FESYSTEM_H
