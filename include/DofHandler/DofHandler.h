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
// define DofHanlder in AsFem

#ifndef ASFEM_DOFHANDLER_H
#define ASFEM_DOFHANDLER_H

#include <iostream>
#include <vector>

#include "petsc.h"

// AsFem's own header file
#include "Mesh/Mesh.h"
#include "BCSystem/BCInfo.h"

using namespace std;

class DofHandler
{
public:
    DofHandler();

    void SetDofsPerNode(int ndofspernode);
    int GetDofsPerNode() const { return nDofsPerNode;}
    int GetDofsNum() const { return nDofs;}

    bool CreateLocalToGlobalDofMap(Mesh &mesh,int ndofspernode);
    bool SetNodalDofActiveState(Mesh &mesh,BCInfo &bcInfo);


    // TODO:implement gmsh dof map generation
    bool CreateLocalToGlobalDofMapFromGmsh(Mesh &mesh,int ndofspernode);


    int GetBCSideDofsNum(string sidename) const;
    int GetBCSideElmtsNum(string sidename) const { return GetBCSideDofsNum(sidename)/nDofsPerBCElmt;}
    double GetIthNodalJthDofState(const int &i,const int &j) const;

    void GetLocalDofMap(const int e,int &ndofsperelmt,int (&rInd)[270]) const;
    void GetLocalBCDofMap(string sidename,const int e,int &ndofsperbcelmt,int (&Ind)[270]) const;

    void Release();

    void PrintDofMap() const;

private:
    bool HasDofMap=false;
    bool HasBCDofMap=false;
    int nDims,nDofs,nNodes,nElmts;
    int nDofsPerNode,nNodesPerElmts;
    int nDofsPerElmt;

    vector<int> GlobalDofMap;
    vector<double> NodalDofState;
    vector<pair<string,vector<int>>> GlobalBCDofMap;
    vector<int> LeftBCDofs,RightBCDofs;
    vector<int> BottomBCDofs,TopBCDofs;
    vector<int> BackBCDofs,FrontBCDofs;

    int nDofsPerBCElmt,nNodesPerBCElmt;
};


#endif //ASFEM_DOFHANDLER_H
