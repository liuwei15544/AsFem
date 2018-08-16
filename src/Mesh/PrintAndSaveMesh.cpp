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
// Define the basic mesh class

#include "Mesh/Mesh.h"

void Mesh::SaveMeshToVTU(string filename)
{
    PetscMPIInt rank,size;
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    if(rank==0)
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
            for(j=1;j<=GetNodesNumPerElmt();++j)
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
            offset+=GetNodesNumPerElmt();
            out<<setw(8)<<offset<<"\n";
        }
        out<<"</DataArray>\n";
        out<<"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n";
        for(e=1;e<=GetElmtsNum();++e)
        {
            out<<setw(8)<<VTKCellType<<"\n";
        }
        out<<"</DataArray>\n";
        out<<"</Cells>\n";
        out<<"</Piece>\n";
        out<<"</UnstructuredGrid>\n";
        out<<"</VTKFile>"<<endl;
        out.close();
    }
}

//*******************************************
void Mesh::PrintMeshInfo(string str) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print mesh info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Mesh information:                      ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nDims=%2d                             ***\n",nDim);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nNodesPerElmt);

    int i,j;
    i=LeftBCConn.size();j=GetSideBCElmtNum("left");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Left side nodes=%5d, elmts=%5d   ***\n",i,j);

    i=RightBCConn.size();j=GetSideBCElmtNum("right");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Right side nodes=%5d, elmts=%5d  ***\n",i,j);

    if(nDim==2)
    {
        i=BottomBCConn.size();j=GetSideBCElmtNum("bottom");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Bottom side nodes=%5d, elmts=%5d ***\n",i,j);

        i=TopBCConn.size();j=GetSideBCElmtNum("top");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   Top side nodes=%5d, elmts=%5d    ***\n",i,j);
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}

//*****************************
void Mesh::PrintMeshDetailInfo(string str) const
{
    int i,j,e;
    double x,y,z,w;
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print mesh info, mesh hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        CreateMesh should be called before this!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }


    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Mesh detailed information:             ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nDims=%2d                             ***\n",nDim);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodes=%10d                    ***\n",nNodes);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nElmts=%10D                    ***\n",nElmts);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   nNodesPerElmt=%3d                    ***\n",nNodesPerElmt);
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Nodes' coordinates:                    ***\n");
    for(i=1;i<=nNodes;++i)
    {
        w=GetIthNodeJthCoord(i,0);
        x=GetIthNodeJthCoord(i,1);
        y=GetIthNodeJthCoord(i,2);
        z=GetIthNodeJthCoord(i,3);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %6d-th node:",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"x=%12.5e,y=%12.5e,z=%12.5e,w=%6.2f\n",x,y,z,w);
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Elements's connectivity:                ***\n");
    for(e=1;e<=nElmts;e++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** e=%6d:",e);
        for(i=1;i<=nNodesPerElmt;++i)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",GetIthConnJthIndex(e,i));
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Boundary element's connectivity:       ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    string sidename;
    for(unsigned int i=1;i<=BCMeshSet.size();++i)
    {
        sidename=BCMeshSet[i-1].first;

        for(e=1;e<=GetSideBCElmtNum(sidename);e++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %8s side->",sidename.c_str());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"e=%3d:",e);
            for(j=1;j<=nNodesPerBCElmt;j++)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",GetSideBCIthConnJthIndex(sidename,e,j));
            }
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}
