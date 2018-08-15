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
// get mesh information of AsFem

#include "Mesh/Mesh.h"

int Mesh::GetIthConnJthIndex(int e, int j) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh is not generated !!!            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!!      ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d !         ***\n",e,nElmts);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    if(j<1||j>nNodesPerElmt)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%3d is out of nNodesPerElmt=%3d !!!***\n",j,nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 1~%3d !!!                ***\n",nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    return Conn[(e-1)*nNodesPerElmt+j-1];
}
//***********************************
double Mesh::GetIthNodeJthCoord(int i, int j) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh is not generated !!!            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!!      ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>nNodes)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%6d is out of nNodes=%6d!          ***\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }
    if(j<0||j>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out 0~3\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 0->w,1->x,2->y,3->z!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return NodeCoords[4*(i-1)+j];
}
//*************************************
void Mesh::GetLocalCoords(int e, double (&coords)[27][4]) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh is not generated !!!       ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,nElmts);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    int i,j;
    for(i=1;i<=nNodesPerElmt;i++)
    {
        j=GetIthConnJthIndex(e,i);
        coords[i-1][0]=GetIthNodeJthCoord(j,0);
        coords[i-1][1]=GetIthNodeJthCoord(j,1);
        coords[i-1][2]=GetIthNodeJthCoord(j,2);
        coords[i-1][3]=GetIthNodeJthCoord(j,3);
    }
}

