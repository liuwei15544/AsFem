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


class FESystem
{
public:
    FESystem();
    void Init(Mesh &mesh);


    // For Ax=F system
    void FormKR(const int &iState,const double dt,const double t,const double (&ctan)[2],
                Mesh &mesh,DofHandler &dofHandler,BCSystem &bcSystem,
                const Vec &U,const Vec &V,
                Mat &AMATRIX,Vec &RHS);
    // For assemble
    void AssembleLocalToGlobal(Mat &A,Vec &RHS);
    void AssembleLocalRHSToGlobal(Vec &RHS);
    void AssembleLocalKToGlobal(Mat &A);

private:
    int nDims,nNodesPerElmt,nDofsPerElmt,nDofsPerNode;
    int current_elmt_id;
    double localHist[50],localProj[12+1];
    double localK[270][270],localRHS[270]; // local matrix and residual
    double elCoords[27][4],elU[270][2]={0.0};
    PetscInt elConn[27]={0},elDofsConn[270]={0};
    double t,dt;

private:
    VecScatter scatteru,scatterv;
    Vec Useq;// this can contain the ghost node from other processor
    Vec Vseq;
};


#endif //ASFEM_FESYSTEM_H
