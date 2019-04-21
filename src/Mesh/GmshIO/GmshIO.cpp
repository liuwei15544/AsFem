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
// Created by walkandthinker on 08.04.19.
// Define the gmsh io class(read .msh file from gmsh)

#include "Mesh/GmshIO.h"
#include "Utils/StringUtils.h"

GmshIO::GmshIO()
{
    MshFileName="";
    version=2.0;
    HasMshFileName=false;

    Xmin=1.0e15;Xmax=-1.0e15;
    Ymin=1.0e15;Ymax=-1.0e15;
    Zmin=1.0e15;Zmax=-1.0e15;
    nDims=0;nElmts=0;nNodes=0;
    nPhysics=0;
    MaxPhyDim=-10;MinPhyDim=10;ElmtMaxDim=0;
    nVolumeElmts=0;nSurfaceElmts=0;nLineElmts=0;nPointElmts=0;
    nBulkElmts=0;

    BulkElmtTypeName="";
}

GmshIO::GmshIO(string mshfilename)
{
    MshFileName=RemoveSpace(mshfilename);
    HasMshFileName=true;
    version=2.0;

    Xmin=1.0e15;Xmax=-1.0e15;
    Ymin=1.0e15;Ymax=-1.0e15;
    Zmin=1.0e15;Zmax=-1.0e15;
    nDims=0;nElmts=0;nNodes=0;
    nPhysics=0;
    MaxPhyDim=-10;MinPhyDim=10;ElmtMaxDim=0;

    nBulkElmts=0;
    nVolumeElmts=0;nSurfaceElmts=0;nLineElmts=0;nPointElmts=0;
    BulkElmtTypeName="";
}

