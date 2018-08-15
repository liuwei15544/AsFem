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

#ifndef SFEM_MESH_H
#define SFEM_MESH_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "petsc.h"

using namespace std;


class Mesh
{
public:
    Mesh();
    void Release();

    bool IsMeshCreated() const { return MeshCreated;}

    // For mesh settings
    void SetDim(int dim);
    void SetXmin(double xmin=0.0);
    void SetXmax(double xmax=1.0);
    void SetYmin(double ymin=0.0);
    void SetYmax(double ymax=1.0);
    void SetZmin(double zmin=0.0);
    void SetZmax(double zmax=1.0);
    void SetMeshType(string meshtype);
    void SetNx(int nx);
    void SetNy(int ny);
    void SetNz(int nz);
    bool IsMeshInfoComplete();

    // For getting mesh information
    int GetDims() const { return nDim;}
    int GetNodesNumPerElmt() const { return nNodesPerElmt;}
    int GetNodesNumPerBCElmt() const { return nNodesPerBCElmt;}
    int GetNodesNum() const { return nNodes;}
    int GetElmtsNum() const { return nElmts;}
    int GetBCElmtsNum() const { return nBCElmts;}


    int GetIthConnJthIndex(int e, int j) const;
    double GetIthNodeJthCoord(int i, int j) const;
    void GetLocalCoords(int e,double (&coords)[27][4]) const;

    //* For boundary mesh information
    int GetSideBCElmtNum(string sidename) const;
    int GetSideBCIthConnJthIndex(string sidename,int i,int j) const;
    double GetSideBCIthNodeJthCoord(string sidename,int e,int i,int j) const;
    void GetLocalBCCoords(string sidename,int e,double (&coords)[27][4]) const;
    //int GetBCIthConnJthIndex(int i,int j);
    //double GetBCIthNodeJthCoord(int i,int j) const;

    //* For VTK information
    //void SetVTKCellType(int vtkcelltype);
    int GetVTKCellType() const { return VTKCellType;}

    // For bult-in mesh
    void CreateMesh();
    void Create1DMesh();
    void Create2DMesh();
    void Create3DMesh();

    // For external mesh
    void ImportGmsh();

    void SaveMeshToVTU(string filename="mesh.vtu");

    void PrintMeshInfo(string str="") const;
    void PrintMeshDetailInfo(string str="") const;


private:
    // for basic mesh information
    double Xmax,Xmin,Ymax,Ymin,Zmax,Zmin;
    bool MeshCreated=false;
    int VTKCellType;
    string MeshType;
    string GmshFileName;
    int Nx,Ny,Nz,nDim;
    int nNodes,nElmts,nNodesPerElmt;
    vector<double> NodeCoords;
    vector<int> Conn;

    // for boundary mesh
    bool BCMeshCreated=false;
    int nBCNodes,nBCElmts,nNodesPerBCElmt;
    vector<int> BCConn;
    vector<int> LeftBCConn,RightBCConn;
    vector<int> BottomBCConn,TopBCConn;
    vector<int> BackBCConn,FrontBCConn;
    vector<pair<string,vector<int>>> BCMeshSet;

private:
    // for state variable check
    bool IsBultInMesh;
    bool IsXminSet,IsXmaxSet;
    bool IsYminSet,IsYmaxSet;
    bool IsZminSet,IsZmaxSet;
    bool IsMeshTypeSet;
    bool IsNxSet,IsNySet,IsNzSet;
    bool IsDimSet;

};


#endif //SFEM_MESH_H
