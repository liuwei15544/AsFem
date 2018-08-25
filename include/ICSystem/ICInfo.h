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
// define the initial condition info for input file(combine all the ic block)

#ifndef ASFEM_ICINFO_H
#define ASFEM_ICINFO_H

#include <iostream>
#include <vector>

//*********************************
#include "ICBlockInfo.h"
#include "EquationSystem/EquationSystem.h"

using namespace std;

class ICInfo
{
public:
    ICInfo();

    void Init();

    int GetICBlockNum() const { return ICBlockList.size();}
    int GetIthICKernelDofIndex(int i) const;
    string GetIthICKernelName(int i) const;
    //double GetIthICKernelValue(int i) const;

    ICBlockInfo GetIthICKernelBlock(int i) const;

    void AddICKernelBlock(ICBlockInfo &icBlockInfo);
    void AddICKernelBlockFromList(vector<ICBlockInfo> icBlockList);

    bool CheckDofsName(EquationSystem &equationSystem);
    void GenerateICKernelDofMap(EquationSystem &equationSystem);

    void PrintICInfo() const;

private:
    bool HasICKernelBlocks=false;
    bool DofsNameCheckPassed=false;
    vector<ICBlockInfo> ICBlockList;
    vector<int> ICBlockDofIndex;
};


#endif //ASFEM_ICINFO_H
