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
// define the boundary block info for input file

#ifndef ASFEM_BCBLOCKINFO_H
#define ASFEM_BCBLOCKINFO_H


#include <iostream>
#include <string>

#include "petsc.h"


using namespace std;

class BCBlockInfo
{
public:
    string BCBlockName;      //[leftbc]-->BCBlockName
    string BCBlockKernelName;// type=dirichlet-->BC kernel name
    string BCBlockDofsName;  // dofs=c        -->BC dofs name
    int    BCBlockDofsIndex; // index of c
    string sidename;         // boundary=left -->BC side name
    double BCValue=0.0;

    void Clear()
    {
        BCBlockName="";
        BCBlockKernelName="";
        BCBlockDofsName="";
        BCBlockDofsIndex=0;
        sidename="";
        BCValue=0.0;
    }

};


#endif //ASFEM_BCBLOCKINFO_H
