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
// print the information for mesh class

#include "Mesh/Mesh.h"


void Mesh::PrintMeshInfo() const
{
    cout<<"*********************************************************"<<endl;
    cout<<"*** Mesh information summary:                         ***"<<endl;
    if(IsBultInMesh)
    {
        cout<<"***   use AsFem's bult-in mesh generation             ***"<<endl;
    }
    else
    {
        if(UseMeshFromGmsh)
        {
            printf("***   use mesh from gmsh(file= %-20s)  ***\n",GmshFileName.c_str());
            printf("***   bulk element type = %-10s                  ***\n",BulkElmtTypeName.c_str());
        }
        if(UseMeshFromAbaqus)
        {
            printf("***   use mesh from abaqus(file= %-20s)***\n",AbaqusFileName.c_str());
        }
    }
    cout<<"***   total elements number =  "<<setw(10)<<GetElmtsNum()
        <<"             ***"<<endl;
    cout<<"***   bulk elements number =   "<<setw(10)<<GetBulkElmtsNum()
        <<"             ***"<<endl;

    if(GetLineElmtsNum()!=0)
    {

        cout<<"***   line elements number =   "<<setw(10)<<GetLineElmtsNum()
            <<"             ***"<<endl;

    }
    if(GetSurfaceElmtsNum()!=0)
    {
        cout<<"***   surface elements number =   "<<setw(10)<<GetSurfaceElmtsNum()
            <<"          ***"<<endl;


    }
    if(GetVolumeElmtsNum()!=0)
    {
        cout<<"***   volume elements number =   "<<setw(10)<<GetVolumeElmtsNum()
            <<"           ***"<<endl;
    }

    cout<<"***   nodes number="<<setw(10)<<GetNodesNum()
        <<"                         ***"<<endl;
    cout<<"***   dimension="<<setw(2)<<nDim;
    cout<<", min elmt dim="<<setw(2)<<MinPhyDim;
    cout<<", max elmt dim="<<setw(2)<<MaxPhyDim<<"  ***"<<endl;
    cout<<"***   physical group number = "<<setw(8)<<PhyGroup.size()
        <<"                ***"<<endl;
    for(unsigned int i=0;i<PhyGroup.size();i++)
    {
        cout<<"***   phy id="<<setw(5)<<PhyIDIndex[i];
        cout<<" ,phy dim="<<setw(2)<<PhyGroup[PhyIDIndex[i]-1].first;
        printf(" ,phy name= %-12s***\n",PhyGroup[PhyIDIndex[i]-1].second.c_str());
    }
    cout<<"*********************************************************"<<endl;
}
