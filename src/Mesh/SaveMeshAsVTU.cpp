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
// Created by walkandthinker on 18.04.19.
// save mesh as vtu file


#include "Mesh/Mesh.h"


void Mesh::SaveMeshAsVTU(string filename)
{
    ofstream out;
    string VTUFileName;
    int i,j,e;

    if(filename.size()>4)
    {
        i=filename.size()-4;
        if(filename.compare(i,4,".vtu")==0)
        {
            VTUFileName=filename;
        }
        else
        {
            VTUFileName=filename+".vtu";
        }
    }
    else
    {
        VTUFileName=filename+".vtu";
    }


    out.open(VTUFileName.c_str(),ios::out);

    out<<"<?xml version=\"1.0\"?>\n";
    out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    out<<"<UnstructuredGrid>\n";
    out<<"<Piece NumberOfPoints=\""<<GetNodesNum()<<"\" NumberOfCells=\""<<GetElmtsNum()<<"\">\n";
    out<<"<Points>\n";
    out<<"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n";

    for(i=1;i<=GetNodesNum();++i)
    {
        out<<scientific<<setprecision(6)
           <<GetIthNodeJthCoord(i,1)<<"  "
           <<GetIthNodeJthCoord(i,2)<<"  "
           <<GetIthNodeJthCoord(i,3)<<"\n";
    }


    out<<"</DataArray>\n";
    out<<"</Points>\n";
    out<<"<Cells>\n";
    out<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(e=1;e<=GetElmtsNum();++e)
    {
        for(j=1;j<=GetIthElmtNodesNum(e);++j)
        {
            out<<" "<<setw(8)<<GetIthConnJthIndex(e,j)-1;
        }
        out<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    long int offset=0;
    for(e=1;e<=GetElmtsNum();++e)
    {
        offset+=GetIthElmtNodesNum(e);
        out<<setw(8)<<offset<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
    for(e=1;e<=GetElmtsNum();++e)
    {
        out<<setw(8)<<GetIthElmtVTKCellType(e)<<"\n";
    }
    out<<"</DataArray>\n";
    out<<"</Cells>\n";
    out<<"</Piece>\n";
    out<<"</UnstructuredGrid>\n";
    out<<"</VTKFile>"<<endl;
    out.close();

}

