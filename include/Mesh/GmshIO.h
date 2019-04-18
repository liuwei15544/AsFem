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

    bool ReadMshFile(vector<double> &NodeCoords,
                     vector<vector<int>> &Conn,
                     vector<pair<int,string>> &GmshPhyGroup);

    void SetMshFileName(string mshfilename) {MshFileName=mshfilename;HasMshFileName=true;}
    void SetMshFileVersion(double ver) {version=ver;}

    int GetNodesNumViaElmtType(int elmttype) const;
    int GetElmtDimViaElmtType(int elmttype) const;
    string GetElmtNameViaElmtType(int elmttype) const;
    vector<int> GetNodeOrderViaElmtType(int elmttype) const;
    int GetVTKCellTypeViaElmtType(int elmttype) const;

    inline int GetMinPhyDim() const { return MinPhyDim;}
    inline int GetMaxPhyDim() const { return MaxPhyDim;}
    inline int GetMaxElmtDim() const { return ElmtMaxDim;}

private:
    bool HasMshFileName=false;
    string MshFileName="";
    double version;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    int nDims,nElmts,nNodes;
    int nPhysics=0;
    int MaxPhyDim=-10,MinPhyDim=10,ElmtMaxDim=0;

};


#endif //ASFEM_GMSHIO_H
