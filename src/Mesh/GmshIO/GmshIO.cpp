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
}

GmshIO::GmshIO(string mshfilename)
{
    MshFileName=RemoveSpace(mshfilename);
}

