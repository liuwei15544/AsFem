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
#include "DofHandler/DofHandler.h"
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

    void ApplyDirichletBC(Mesh &mesh,DofHandler &dofHandler,
                          const double &t,const double &dt,Vec &U);

    void ApplyNeumannBC(Mesh &mesh,DofHandler &dofHandler,
                        const double &t,const double &dt,
                        Mat &AMATRIX,Vec &RHS);
    void ubc1(string sidename,const double value,const int &dofindex,
              Mesh &mesh,DofHandler &dofHandler,
              const double &t,const double &dt,
              Mat &AMATRIX,Vec &RHS);

    void ApplyConstraint(Mesh &mesh,DofHandler &dofHandler,Mat &AMATRIX,Vec &RHS);

    void PrintBCInfo() const;
    BCInfo bcInfo;

private:
    BCBlockInfo SingleBCBlock;

    bool IsInit=false;
    int nDims;

    PetscMPIInt rank,size;
};


#endif //ASFEM_BCSYSTEM_H
