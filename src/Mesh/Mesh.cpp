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

#include "Mesh/Mesh.h"

Mesh::Mesh()
{
    // for basic mesh information

    Xmax=0.0;Xmin=1.0;
    Ymax=0.0;Ymin=1.0;
    Zmax=0.0;Zmin=1.0;

    MeshCreated=false;
    MeshType="";
    GmshFileName="";
    Nx=0;Ny=0;Nz=0;nDim=0;
    nNodes=0;nElmts=0;nNodesPerElmt=0;
    NodeCoords.clear();
    Conn.clear();

    // for boundary mesh
    BCMeshCreated=false;
    nBCNodes=0;nBCElmts=0;nNodesPerBCElmt=0;

    BCMeshSet.clear();
    BCConn.clear();
    LeftBCConn.clear();RightBCConn.clear();
    BottomBCConn.clear();TopBCConn.clear();
    BackBCConn.clear();FrontBCConn.clear();

    // for state variable check
    IsBultInMesh=true;
    IsXminSet=false;IsXmaxSet=false;
    IsYminSet=false;IsYmaxSet=false;
    IsZminSet=false;IsZmaxSet=false;
    IsMeshTypeSet=false;
    IsNxSet=false;IsNySet=false;IsNzSet=false;
    IsDimSet=false;
}

//*************************
void Mesh::Release()
{
    NodeCoords.clear();
    Conn.clear();

    BCMeshSet.clear();
    BCConn.clear();
    LeftBCConn.clear();RightBCConn.clear();
    BottomBCConn.clear();TopBCConn.clear();
    BackBCConn.clear();FrontBCConn.clear();
}

//***********************************
bool Mesh::IsMeshInfoComplete()
{
    if(!IsDimSet)
    {
        return false;
    }
    else
    {
        if(nDim==1)
        {
            if(!IsNxSet) return false;
            if(!IsXminSet) return false;
            if(!IsXmaxSet) return false;
        }
        else if(nDim==2)
        {
            if(!IsNxSet) return false;
            if(!IsNySet) return false;

            if(!IsXminSet) return false;
            if(!IsXmaxSet) return false;
            if(!IsYminSet) return false;
            if(!IsYmaxSet) return false;
        }
        else if(nDim==3)
        {
            if(!IsNxSet) return false;
            if(!IsNySet) return false;
            if(!IsNzSet) return false;

            if(!IsXminSet) return false;
            if(!IsXmaxSet) return false;
            if(!IsYminSet) return false;
            if(!IsYmaxSet) return false;
            if(!IsZminSet) return false;
            if(!IsZmaxSet) return false;
        }
    }

    if(!IsMeshTypeSet) return false;

    return true;
}

