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
#include <map>
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
                     vector<int> &ElmtVTKCellType,
                     vector<string> &ElmtTypeName,
                     map<string,vector<int>> &MeshNameSet,
                     map<int,vector<int>> &MeshIdSet,
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
    inline int GetElmtsNum() const { return nElmts;}
    inline int GetBulkElmtsNum() const { return nBulkElmts;}
    inline int GetPointElmtsNum() const { return nPointElmts;}
    inline int GetLineElmtsNum() const { return nLineElmts;}
    inline int GetSurfaceElmtsNum() const { return nSurfaceElmts;}
    inline int GetVolumeElmtsNum() const { return nVolumeElmts;}
    inline int GetNodesNum() const { return nNodes;}

    inline double GetXmin() const { return Xmin;}
    inline double GetXmax() const { return Xmax;}

    inline double GetYmin() const { return Ymin;}
    inline double GetYmax() const { return Ymax;}

    inline double GetZmin() const { return Zmin;}
    inline double GetZmax() const { return Zmax;}

private:
    bool HasMshFileName=false;
    string MshFileName="";
    double version;
    double Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;
    int nDims=0,nElmts=0,nNodes=0;
    int nPhysics=0;
    int MaxPhyDim=-10,MinPhyDim=10,ElmtMaxDim=0;
    int nVolumeElmts,nSurfaceElmts,nLineElmts,nPointElmts;
    int nBulkElmts;



};


#endif //ASFEM_GMSHIO_H
