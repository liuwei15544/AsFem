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
// Generate 1D mesh for AsFem

#include "Mesh/Mesh.h"

void Mesh::Create1DMesh()
{
    MeshCreated=false;
    int i,j,e,P;

    VTKCellType=4;P=1;
    if(MeshType=="edge2")
    {
        P=1;
        VTKCellType=3;
    }
    if(MeshType=="edge3") P=2;
    if(MeshType=="edge4") P=3;

    nElmts=Nx;
    nNodesPerElmt=P+1;
    nNodes=nElmts*P+1;

    double dx=(Xmax-Xmin)/(nNodes-1.0);


    NodeCoords.resize(nNodes*4,0.0);
    Conn.resize(nElmts*nNodesPerElmt,0);




    for(i=0;i<nNodes;++i)
    {
        NodeCoords[i*4  ]=1.0;
        NodeCoords[i*4+1]=Xmin+i*dx;
        NodeCoords[i*4+2]=0.0;
        NodeCoords[i*4+3]=0.0;
    }


    for(e=0;e<nElmts;++e)
    {
        for(j=1;j<=nNodesPerElmt;++j)
        {
            Conn[e*nNodesPerElmt+j-1]=e*P+j;
        }
    }

    //***********************************
    //*** generate boundary mesh info
    //***********************************
    BCConn.resize(2,0);
    nBCElmts=2;nNodesPerBCElmt=1;

    BCConn[0]=1;BCConn[1]=nNodes;
    LeftBCConn.resize(1);LeftBCConn[0]=1;
    RightBCConn.resize(1);RightBCConn[0]=nNodes;
    BottomBCConn.clear();TopBCConn.clear();
    BackBCConn.clear();FrontBCConn.clear();
    BCMeshSet.clear();

    pair<string,vector<int>> temppair;

    temppair=make_pair("left",LeftBCConn);
    BCMeshSet.push_back(temppair);

    temppair=make_pair("right",RightBCConn);
    BCMeshSet.push_back(temppair);

    MeshCreated=true;
    BCMeshCreated=true;
}

