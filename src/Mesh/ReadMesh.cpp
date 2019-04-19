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
// Created by walkandthinker on 19.04.19.
// Read mesh from mesh file


#include "Mesh/Mesh.h"

bool Mesh::ReadMesh(string meshfiletype)
{
    if(meshfiletype=="gmsh")
    {
        GmshIO gmshIo(GmshFileName);

        if(gmshIo.ReadMshFile(NodeCoords,
                           Conn,
                           MeshNameSet,
                           MeshIdSet,
                           PhyGroup))
        {
            Xmin=gmshIo.GetXmin();
            Xmax=gmshIo.GetXmax();

            Ymin=gmshIo.GetYmin();
            Ymax=gmshIo.GetYmax();

            Zmin=gmshIo.GetZmin();
            Zmax=gmshIo.GetZmax();
            return true;
        }
        else
        {
            Msg_Gmsh_ImportFailed();
            return false;
        }
    }
    else if(meshfiletype=="abaqus")
    {
        cout<<"*** Sorry, abaqus mesh import is not supported yet!"<<endl;
        return false;
    }
}

