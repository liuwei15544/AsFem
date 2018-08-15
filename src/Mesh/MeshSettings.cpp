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
// Define the settings for mesh class

#include "Mesh/Mesh.h"

void Mesh::SetDim(int dim)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid for FEM mesh!!!   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    nDim=dim;
    IsDimSet=true;
}

//***************************************
void Mesh::SetXmin(double xmin)
{
    if(IsXmaxSet)
    {
        if(xmin>=Xmax)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: xmin=%9.3f >= xmax=%9.3f !!!   ***\n",xmin,Xmax);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    Xmin=xmin;
    IsXminSet=true;
}
void Mesh::SetXmax(double xmax)
{
    if(IsXminSet)
    {
        if(xmax<=Xmin)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: xmax=%9.3f <= Xmin=%9.3f !!!   ***\n",xmax,Xmin);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    Xmax=xmax;
    IsXmaxSet=true;
}
//********************************
void Mesh::SetYmin(double ymin)
{
    if(IsDimSet)
    {
        if(nDim<2)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cant' set Ymin, dim=%3d !!!            ***\n",nDim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    if(IsYmaxSet)
    {
        if(ymin>=Ymax)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ymin=%9.3f >= Ymin=%9.3f !!!   ***\n",ymin,Ymax);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    Ymin=ymin;
    IsYminSet=true;
}
void Mesh::SetYmax(double ymax)
{
    if(IsDimSet)
    {
        if(nDim<2)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cant' set Ymax, dim=%3d !!!            ***\n",nDim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    if(IsYminSet)
    {
        if(ymax<=Ymin)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ymax=%9.3f <= Ymin=%9.3f !!!   ***\n",ymax,Ymin);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    Ymax=ymax;
    IsYmaxSet=true;
}
//*******************************
void Mesh::SetZmin(double zmin)
{
    if(IsDimSet)
    {
        if(nDim<3)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cant' set Zmin, dim=%3d !!!            ***\n",nDim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(IsZmaxSet)
    {
        if(zmin>=Zmax)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: zmin=%9.3f >= Zmax=%9.3f !!!   ***\n",zmin,Zmax);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }
    Zmin=zmin;
    IsZminSet=true;
}
void Mesh::SetZmax(double zmax)
{
    if(IsDimSet)
    {
        if(nDim<3)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cant' set Zmax, dim=%3d !!!            ***\n",nDim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(IsZminSet)
    {
        if(zmax<=Zmin)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: zmax=%9.3f <= Zmin=%9.3f !!!   ***\n",zmax,Zmin);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }

    Zmax=zmax;
    IsZmaxSet=true;
}

//*************************************
//*** Mesh number setting
//*************************************
void Mesh::SetNx(int nx)
{
    if(nx<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%6d is invalid                   ***\n",nx);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscFinalize();
        abort();
    }
    Nx=nx;
    IsNxSet=true;
}
