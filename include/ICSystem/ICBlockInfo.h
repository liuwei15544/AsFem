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
// define initial condition block info for input file

#ifndef ASFEM_ICBLOCKINFO_H
#define ASFEM_ICBLOCKINFO_H

#include <iostream>
#include <vector>
#include <string>

#include "petsc.h"


using namespace std;

class ICBlockInfo
{
public:
    string ICBlockName;      //[leftbc]-->BCBlockName
    string ICBlockKernelName;// type=dirichlet-->BC kernel name
    string ICBlockDofsName;  // dofs=c        -->BC dofs name
    int    ICBlockDofsIndex; // index of c
    vector<double> ICParams;

    void Clear()
    {
        ICBlockName="";
        ICBlockKernelName="";
        ICBlockDofsName="";
        ICBlockDofsIndex=0;
        ICParams.clear();
    }
};


#endif //ASFEM_ICBLOCKINFO_H
