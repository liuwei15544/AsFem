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
// Created by walkandthinker on 17.11.18.
// import mesh from .msh file(by Gmsh)

#include "Mesh/Mesh.h"

void Mesh::ImportGmsh(string filename)
{
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if(rank==0)
    {
        ifstream in;

        if(filename.length()<5)
        {
            // at least must be '*.msh'
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: msh file name(=%s) is invalid!!!\n",filename.c_str());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
            PetscFinalize();
            abort();
        }

        GmshFileName=filename;
        in.open(GmshFileName.c_str(),ios::in);
        string line,str,str0;

        int phydim,phyid;
        string phyname;


        while (!in.eof())
        {
            getline(in,line);
            if(line.find("$PhysicalNames")!=string::npos)
            {
                nPhysics=0;
                MaxPhyDim=-10;MinPhyDim=10;
                in>>nPhysics;
                GmshPhyGroup.clear();
                GmshPhyGroup.resize(nPhysics);
                for(int i=0;i<nPhysics;i++)
                {
                    //phyical dimension,physical id,physical name
                    getline(in,line);
                    istringstream s_stream(line);
                    s_stream >> phydim >> phyid >> phyname;

                    if(phydim>MaxPhyDim) MaxPhyDim=phydim;
                    if(phydim<MinPhyDim) MinPhyDim=phydim;

                    //remove ""
                    phyname.erase(remove(phyname.begin(),phyname.end(),'"'), phyname.end());
                    GmshPhyGroup[phyid-1]=make_pair(phydim,phyname);
                }

                getline(in,line);// read '$EndPhysicalNames'
            }
            else if(line.find("$NOD")!=string::npos ||
                    line.find("$NOE")!=string::npos||
                    line.find("$Nodes")!=string::npos)
            {
                // Read node's coordinates
                nNodes=0;
                in>>nNodes;
                NodeCoords.clear();
                NodeCoords.resize(4*nNodes,0.0);

                int nodeid;
                double x,y,z;

                for(int i=0;i<nNodes;i++)
                {
                    in>>nodeid>>x>>y>>z;
                    NodeCoords[(nodeid-1)*4+0]=1.0;
                    NodeCoords[(nodeid-1)*4+1]=x;
                    NodeCoords[(nodeid-1)*4+2]=y;
                    NodeCoords[(nodeid-1)*4+3]=z;
                }

                getline(in,line);// read '$EndNodes'
            }
            else if(line.find("$Elements")!=string::npos||
                    line.find("$ELM")!=string::npos)
            {
                in>>nElmts;
                Conn.clear();
                //elmt-id, elmt-type, tag num,
                int elmtid,elmttype,ntags;
                int nBCElmt=0;
                for(int i=0;i<nElmts;i++)
                {
                    in>>elmtid>>elmttype>>ntags;
                    if(GetNodesNumOfGmshCell(elmttype)==1)
                    {
                        nBCElmt+=1;
                    }
                    else
                    {

                    }
                }
            }
        }
    }
}

//**************************************************
int Mesh::GetNodesNumOfGmshCell(int elmttpye)
{
    if(elmttpye==1)
    {
        // 2-node line
        return 2;
    }
    else if(elmttpye==8)
    {
        // 3-node second order line
        // 1--3--2
        return 3;
    }
    else if(elmttpye==26)
    {
        // 4-node third order edge
        // 1--3--4--2
        return 4;
    }
    else if(elmttpye==2)
    {
        // 3-node triangle
        return 3;
    }
    else if(elmttpye==9)
    {
        // 6-node second order triangle
        return 6;
    }
    else if(elmttpye==3)
    {
        // 4-node quadrangle.
        return 4;
    }
    else if(elmttpye==10)
    {
        // 9-node second order quadrangle
        return 9;
    }
    else if(elmttpye==16)
    {
        // 8-node second order quadrangle
        return 8;
    }
}