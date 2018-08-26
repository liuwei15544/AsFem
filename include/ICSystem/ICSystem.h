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
// define the initial condition system in AsFem

#ifndef ASFEM_ICSYSTEM_H
#define ASFEM_ICSYSTEM_H

#include <iostream>
#include <string>

//*********************************
// For AsFem's own header file
//*********************************
#include "ICInfo.h"
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"

class ICSystem
{
public:
    ICSystem();
    void Init();
    void SetDofsNumPerNode(int ndofspernode) {nDofsPerNode=ndofspernode;}

    // For initializing
    void AddICKernelBlock(ICBlockInfo &icBlockInfo);
    void InitFromICBlockList(vector<ICBlockInfo> &icBlockList);

    void SetUpICSystem(EquationSystem &equationSystem);

    void ApplyInitialCondition(Mesh &mesh,Vec &U);



    void PrintBCInfo() const;

private:
    void ConstantIC(ICBlockInfo &icBlock,Mesh &mesh,Vec &U);
    void CircleIC(ICBlockInfo &icBlock,Mesh &mesh,Vec &U);
    void RandomIC(ICBlockInfo &icBlock,Mesh &mesh,Vec &U);
    void RandomNoiseIC(ICBlockInfo &icBlock,Mesh &mesh,Vec &U);

private:
    int nDofsPerNode=1;
    ICBlockInfo SingleICBlock;
    ICInfo icInfo;
    bool IsInit=false;
    int nDims;

    PetscMPIInt rank,size;
    PetscRandom    rnd;
};


#endif //ASFEM_ICSYSTEM_H
