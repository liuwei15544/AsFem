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

#ifndef ASFEM_GMSHIO_H
#define ASFEM_GMSHIO_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

#include "MsgPrint/MsgPrintForProgram.h"
#include "MsgPrint/MsgPrintForInput.h"

class GmshIO
{
public:
    GmshIO();
    GmshIO(string mshfilename);

    void SetMshFileName(string mshfilename) {MshFileName=mshfilename;}
    void SetMshFileVersion(double ver) {version=ver;}

    int GetNodesNumViaElmtType(int elmttype) const;
    int GetElmtDimViaElmtType(int elmttype) const;
    string GetElmtNameViaElmtType(int elmttype) const;
    vector<int> GetNodeOrderViaElmtType(int elmttype) const;
    int GetVTKCellTypeViaElmtType(int elmttype) const;

private:
    string MshFileName="";
    double version;

};


#endif //ASFEM_GMSHIO_H
