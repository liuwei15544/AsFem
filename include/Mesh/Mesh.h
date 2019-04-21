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

#include "Mesh/GmshIO.h"

class Mesh
{
public:
    Mesh();
    void Release();

    bool CreateMesh();
    bool ReadMesh(string meshfiletype="gmsh");

    inline bool IsMeshCreated() const {return MeshCreated;}

    inline int GetElmtVTKCellType(int e) const { return ElmtVTKCellType[e-1];}
    inline string GetElmtTypeName(int e) const { return ElmtTypeName[e-1];}

    inline int GetNodesNum() const { return nNodes;}
    inline int GetElmtsNum() const { return nElmts;}
    inline int GetBulkElmtsNum() const { return nBulkElmts;}
    inline int GetPointElmtsNum() const { return nPointElmtsNum;}
    inline int GetLineElmtsNum() const { return nLineElmtsNum;}
    inline int GetSurfaceElmtsNum() const { return nSurfaceElmtsNum;}
    inline int GetVolumeElmtsNum() const { return nVolumeElmtsNum;}


    void PrintMeshInfo() const;
    void PrintMeshDetailedInfo() const;

private:
    bool Create1DMesh();
    bool Create2DMesh();
    bool Create3DMesh();
    // for mesh's private variable
private:
    // for basic mesh information
    double Xmax,Xmin,Ymax,Ymin,Zmax,Zmin;
    bool MeshCreated=false;
    int VTKCellType;
    string MeshType;
    int Nx,Ny,Nz,nDim;
    int nNodes,nElmts,nBulkElmts,nNodesPerElmt;
    int nPointElmtsNum,nLineElmtsNum,nSurfaceElmtsNum,nVolumeElmtsNum;
    vector<double> NodeCoords;
    vector<vector<int>> Conn;
    vector<int> ElmtVTKCellType;
    vector<string> ElmtTypeName;

    // For gmsh infortion
    string GmshFileName;
    int nPhysics;
    vector<pair<int,string>> PhyGroup;
    int MaxPhyDim=-10,MinPhyDim=10,ElmtMaxDim=0;
    map<string,vector<int>> MeshNameSet;
    map<int,vector<int>> MeshIdSet;
    // For abaqus
    string AbaqusFileName;

    bool UseMeshFromGmsh=false,UseMeshFromAbaqus=false;

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
