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
// define the FEsystem in AsFem

#ifndef ASFEM_FESYSTEM_H
#define ASFEM_FESYSTEM_H

#include <iostream>
#include <iomanip>

#include "petsc.h"

// For AsFem's own header file
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "EquationSystem/EquationSystem.h"
#include "BCSystem/BCSystem.h"
#include "ElementSystem/ElementSystem.h"

using namespace std;


class FE
{
public:
    FE();
    void Init(Mesh &mesh,DofHandler &dofHandler);


    // Init related vector and matrix
    void ZeroMatAndVec(const int &iState,Mat &AMATRIX,Mat &Proj,Vec &RHS);


    // For Ax=F system
    void FormKR(const int &iState,const double dt,const double t,const double (&ctan)[2],
                Mesh &mesh,DofHandler &dofHandler,ElementSystem &elementSystem,
                const Vec &U,const Vec &V,
                Mat &AMATRIX,Vec &RHS,Mat &Proj);

    // For assemble
    void AssembleLocalToGlobal(const int &iState,Mat &AMATRIX,Vec &RHS,Mat &Proj);
    void FinishAssemble(const int &iState,Mat &AMATRIX,Vec &RHS,Mat &Proj);

private:
    bool IsInit;
    int nDims,nNodesPerElmt,nDofsPerElmt,nDofsPerNode;
    int current_elmt_id;
    double localHist[50],localProj[27][12+1];
    double localK[270][270],localRHS[270]; // local matrix and residual
    double elCoords[27][4],elU[270][2]={0.0};
    PetscInt elConn[27]={0},elDofsConn[270]={0};

    PetscMPIInt rank,size;


private:
    VecScatter scatteru,scatterv;
    Vec Useq;// this can contain the ghost node from other processor
    Vec Vseq;
};


#endif //ASFEM_FESYSTEM_H
