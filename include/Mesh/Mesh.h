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
    void SetMshFileName(string meshfilename) {GmshFileName=meshfilename;}
    void SetInpFileName(string meshfilename) {AbaqusFileName=meshfilename;}

    inline bool IsMeshCreated() const {return MeshCreated;}

    inline int GetIthElmtVTKCellType(int e) const { return ElmtVTKCellType[e-1];}
    inline string GetIthElmtTypeName(int e) const { return ElmtTypeName[e-1];}

    inline int GetDim() const { return nDim;}
    inline int GetNodesNum() const { return nNodes;}
    inline int GetElmtsNum() const { return nElmts;}
    inline int GetBulkElmtsNum() const { return nBulkElmts;}
    inline int GetPointElmtsNum() const { return nPointElmtsNum;}
    inline int GetLineElmtsNum() const { return nLineElmtsNum;}
    inline int GetSurfaceElmtsNum() const { return nSurfaceElmtsNum;}
    inline int GetVolumeElmtsNum() const { return nVolumeElmtsNum;}

    inline int GetXmax() const { return Xmax;}
    inline int GetXmin() const { return Xmin;}
    inline int GetYmax() const { return Ymax;}
    inline int GetYmin() const { return Ymin;}
    inline int GetZmax() const { return Zmax;}
    inline int GetZmin() const { return Zmin;}

    inline int GetElmtsNumViaPhyID(int phyid) const
    {
        for(auto it=MeshIdSet.begin();it!=MeshIdSet.end();it++)
        {
            if(it->first==phyid) return it->second.size();
        }
    }
    inline int GetElmtsNumViaName(string phyname) const
    {
        for(auto it=MeshNameSet.begin();it!=MeshNameSet.end();it++)
        {
            if(it->first==phyname) return it->second.size();
        }
    }

    // For local mesh information
    inline double GetIthNodeJthCoord(int i,int j) const { return NodeCoords[(i-1)*4+j];}
    inline int GetIthConnJthIndex(int e,int j) const { return Conn[e-1][j];}
    inline int GetIthElmtNodesNum(int e) const { return Conn[e-1][0];}
    inline void GetIthElmtConn(const int &e,int &nnodes,int (&elConn)[27])
    {
        nnodes=Conn[e-1][0];
        for(int i=0;i<nnodes;i++)
        {
            elConn[i]=Conn[e-1][1+i];
        }
    }
    inline void GetIthElmtConnViaID(const int &e,const int &phyid,int &nnodes,int (&elConn)[27])
    {
        int i=MeshIdSet[phyid][e-1];
        nnodes=Conn[i-1][0];
        for(int j=0;j<nnodes;j++)
        {
            elConn[j]=Conn[i-1][1+j];
        }
    }
    inline void GetIthElmtConnViaName(const int &e,string phyname,int &nnodes,int (&elConn)[27])
    {
        int i=MeshNameSet[phyname][e-1];
        nnodes=Conn[i-1][0];
        for(int j=0;j<nnodes;j++)
        {
            elConn[j]=Conn[i-1][1+j];
        }
    }



    void PrintMeshInfo() const;
    void PrintMeshDetailedInfo() const;

    void SaveMeshAsVTU(string filename="mesh.vtu");

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
    vector<vector<int>> Conn; // Conn[e][0]--> number of elemental nodes
                              // Conn[e][i]--> i-th node index of e-th element
    vector<int> ElmtVTKCellType;
    vector<string> ElmtTypeName;
    vector<int> PhyIDIndex;// physical id could be disordered
                           // i.e. PhyIDIndex[0]=55---> 1st phyid is 55
                           //      PhyIDIndex[1]=67---> 2nd phyid is 67

    // For gmsh infortion
    string GmshFileName;
    int nPhysics;
    vector<pair<int,string>> PhyGroup;
    int MaxPhyDim=-10,MinPhyDim=10,ElmtMaxDim=0;
    map<string,vector<int>> MeshNameSet;
    map<int,vector<int>> MeshIdSet;
    // For abaqus
    string AbaqusFileName;

    string BulkElmtTypeName="";

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
