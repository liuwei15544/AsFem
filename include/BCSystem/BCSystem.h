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
// define the boundary system in AsFem

#ifndef ASFEM_BCSYSTEM_H
#define ASFEM_BCSYSTEM_H

#include <iostream>
#include <string>

//*********************************
// For AsFem's own header file
//*********************************
#include "BCInfo.h"
#include "Mesh/Mesh.h"
//#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"

class BCSystem
{
public:
    BCSystem();
    void Init();

    // For initializing
    void SetDims(int dim);
    void AddBCKernelBlock(BCBlockInfo &bcBlockInfo);
    void InitFromBCBlockList(vector<BCBlockInfo> &bcBlockList);

    void SetUpBCSystem(EquationSystem &equationSystem);

    void ApplyBoundaryCondition(Mesh &mesh,DofHandler &dofHandler,EquationSystem &equationSystem);
    void ApplyDirichletBC(string sidename,int dofindex,double value,Mesh &mesh,DofHandler &dofHandler,EquationSystem &equationSystem);
    void ApplyNeumannBC(Mesh &mesh,DofHandler &dofHandler,Vec &RHS);

    void ApplyConstraint(Mesh &mesh,DofHandler &dofHandler,EquationSystem &equationSystem);



    void PrintBCInfo() const;
private:
    BCBlockInfo SingleBCBlock;
    BCInfo bcInfo;
    bool IsInit=false;
    int nDims;
};


#endif //ASFEM_BCSYSTEM_H
