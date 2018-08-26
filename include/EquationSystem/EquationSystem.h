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
// Define equation system in AsFem

#ifndef ASFEM_EQUATIONSYSTEM_H
#define ASFEM_EQUATIONSYSTEM_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "petsc.h"

using namespace std;

class EquationSystem
{
public:
    EquationSystem(const int dofs,int dofspernode=1);
    EquationSystem();

    PetscErrorCode Init();

    // get information
    int GetNodalSolutionNum() const { return nDofsPerNode;}
    int GetDofsNumPerNode() const { return nDofsPerNode;}
    int GetDofsNum() const{ return nDofs;}

    string GetDofNameByIndex(const int i) const;
    int GetDofIndexByName(string dofname) const;

    // set solution name and index
    void AddSolutionNameAndIndex(string name,int order);
    void SetSolutionNameFromVector(vector<string> names,vector<int> orders);
    void SetSolutionName();
    void SetDofsNum(int ndofs) {nDofs=ndofs;IsInit=false;}
    void SetDofsNumPerNode(int ndofspernode) {nDofsPerNode=ndofspernode;}


    // update and reinitializing functions
    void UpdateUplusdU();
    void UpdateU0(Vec &u);
    void UpdateV(Vec &v);
    void ReInitKandR();
    void ReInitVec(Vec &v);
    void ReInitMat(Mat &a);

    // for petsc based math utils
    double GetRNorm();
    double GetdUNorm();
    double GetUNorm();
    double GetEnergyNorm();
    double GetVecNorm(Vec &v);
    double GetKNrom();
    double GetMatNorm(Mat &a);

    // release memory
    PetscErrorCode Release();

    void PrintSolutionNameMap(string str="") const;


public:
    Mat AMATRIX;
    Vec RHS,dU,U,V;
    Vec U0;
    Vec Proj;

private:
    const int MaxDofs=3000000;//I want to limite maximum dofs under 300W
    bool IsInit=false;
    PetscErrorCode ierr;
    int nDofs,nDofsPerNode;
    vector<string> solution_name_list;
    vector<int> solution_index_list;
    vector<pair<string,int>> solution_name_map;
    bool SolutionHasName=false;
};


#endif //ASFEM_EQUATIONSYSTEM_H
