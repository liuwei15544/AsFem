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
// define the boundary info for input file(combine all the bc block)

#ifndef ASFEM_BCINFO_H
#define ASFEM_BCINFO_H

#include <iostream>
#include <vector>

//*********************************
#include "BCBlockInfo.h"
#include "EquationSystem/EquationSystem.h"

using namespace std;


class BCInfo
{
public:
    BCInfo();

    void Init();

    int GetBCBlockNum() const { return BCBlockList.size();}
    int GetIthBCKernelDofIndex(int i) const;
    string GetIthBCKernelName(int i) const;
    double GetIthBCKernelValue(int i) const;
    string GetIthBCKernelSideName(int i) const;

    BCBlockInfo GetIthBCKernelBlock(int i) const;

    void AddBCKernelBlock(BCBlockInfo &bcBlockInfo);
    void AddBCKernelBlockFromList(vector<BCBlockInfo> bcBlockList);

    bool CheckDofsName(EquationSystem &equationSystem);
    void GenerateBCKernelDofMap(EquationSystem &equationSystem);

    void PrintBCInfo() const;

private:
    bool HasBCKernelBlocks=false;
    bool DofsNameCheckPassed=false;
    vector<BCBlockInfo> BCBlockList;
    vector<int> BCBlockDofIndex;
};


#endif //ASFEM_BCINFO_H
