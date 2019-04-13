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
// Created by walkandthinker on 15.08.18.
// Define the basic mesh class

#ifndef ASFEM_MESH_H
#define ASFEM_MESH_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

#include "MsgPrint/MsgPrintForProgram.h"

class Mesh
{
public:
    Mesh();
    void Release();

    inline bool IsMeshCreated() const {return MeshCreated;}


    // for mesh's private variable
private:
    // for basic mesh information
    double Xmax,Xmin,Ymax,Ymin,Zmax,Zmin;
    bool MeshCreated=false;
    int VTKCellType;
    string MeshType;
    int Nx,Ny,Nz,nDim;
    int nNodes,nElmts,nNodesPerElmt;
    vector<double> NodeCoords;
    vector<int> Conn;

    // For gmsh infortion
    string GmshFileName;
    int nPhysics;
    vector<pair<int,string>> GmshPhyGroup;
    int MaxPhyDim=-10,MinPhyDim=10;

    // For state variable check
    bool IsBultInMesh=false;
    bool IsXminSet=false,IsXmaxSet=false;
    bool IsYminSet=false,IsYmaxSet=false;
    bool IsZminSet=false,IsZmaxSet=false;
    bool IsMeshTypeSet=false;
    bool IsNxSet=false,IsNySet=false,IsNzSet=false;
    bool IsDimSet=false;
};


#endif //ASFEM_MESH_H
