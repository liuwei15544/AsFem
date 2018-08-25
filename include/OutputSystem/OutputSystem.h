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
// Created by walkandthinker on 25.08.18.
// Define the output system in AsFem

#ifndef ASFEM_OUTPUTSYSTEM_H
#define ASFEM_OUTPUTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"


//***********************************
//*** For AsFem's own header file
//***********************************
#include "Mesh/Mesh.h"
#include "EquationSystem/EquationSystem.h"

using namespace std;

class OutputSystem
{
public:
    OutputSystem();
    void SetInputFileName(string filename);
    string GetInputFileName() const { return InputFileName;}

    void WriteUToVTUFile(Mesh &mesh,EquationSystem &equationSystem,Vec &U);
    void WriteUToVTUFile(Mesh &mesh,EquationSystem &equationSystem,Vec &U,int step);

    void WriteUAndProjToVTUFile(Mesh &mesh,EquationSystem &equationSystem,Vec &U,Vec &Proj);

private:
    bool IsInit;
    string InputFileName;

private:
    VecScatter scatter,scatterproj;
    Vec Useq,PROJseq;
    FILE *fd;
    PetscMPIInt rank;

};


#endif //ASFEM_OUTPUTSYSTEM_H
