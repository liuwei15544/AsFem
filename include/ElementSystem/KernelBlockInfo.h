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
// define the kernel block for input file reader


#ifndef ASFEM_KERNELBLOCKINFO_H
#define ASFEM_KERNELBLOCKINFO_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "petsc.h"

using namespace std;

class KernelBlockInfo
{
public:
    KernelBlockInfo();
    string ElementName;
    string MaterialKernelName;
    vector<double> params;// directly read from input file
    string strain="small";

    //*********************************
    //*** functions
    void Init();
    int GetParamNum() const { return params.size();}
    void PrintKernelBlockInfo() const;
};


#endif //ASFEM_KERNELBLOCKINFO_H
