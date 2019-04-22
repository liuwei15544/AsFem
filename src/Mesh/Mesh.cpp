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
// Created by walkandthinker on 21.04.19.
// define the mesh class for asfem


#include "Mesh/Mesh.h"

Mesh::Mesh()
{
    Xmax=-1.0e8;Xmin=1.0e8;
    Ymax=-1.0e8;Ymin=1.0e8;
    Zmax=-1.0e8;Zmin=1.0e8;

    MeshCreated=false;
    VTKCellType=-1;
    MeshType="";
    Nx=0;Ny=0;Nz=0;
    nDim=0;
    nNodes=0;nElmts=0;nBulkElmts=0;nNodesPerElmt=0;
    nPointElmtsNum=0;nLineElmtsNum=0;
    nSurfaceElmtsNum=0;nVolumeElmtsNum=0;
    NodeCoords.clear();
    Conn.clear();
    ElmtVTKCellType.clear();
    ElmtTypeName.clear();

    PhyIDIndex.clear();

    // For gmsh infortion
    GmshFileName="";
    nPhysics=-1;
    PhyGroup.clear();
    MaxPhyDim=-10;MinPhyDim=10;ElmtMaxDim=0;
    MeshNameSet.clear();
    MeshIdSet.clear();
    // For abaqus
    AbaqusFileName="";

    BulkElmtTypeName="";

    UseMeshFromGmsh=false;UseMeshFromAbaqus=false;

    // For state variable check
    IsBultInMesh=false;
    IsXminSet=false;IsXmaxSet=false;
    IsYminSet=false;IsYmaxSet=false;
    IsZminSet=false;IsZmaxSet=false;
    IsMeshTypeSet=false;
    IsNxSet=false;IsNySet=false;IsNzSet=false;
    IsDimSet=false;
}

