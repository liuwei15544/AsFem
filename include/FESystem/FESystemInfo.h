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
// Created by walkandthinker on 24.08.18.
// define the FESystem info in AsFem

#ifndef ASFEM_FESYSTEMINFO_H
#define ASFEM_FESYSTEMINFO_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "petsc.h"

using namespace std;

class FESystemInfo
{
public:
    FESystemInfo();

    void Init();
    void PrintFESystemInfo() const;

    double current_dt,old_dt;
    double current_time;
    int totalstep;
    int iState;

    bool IsProjOutput;

    string jobtype;
    string inputfilename;

};


#endif //ASFEM_FESYSTEMINFO_H
