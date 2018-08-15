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
//*********************
void Mesh::SetNy(int ny)
{
    if(IsDimSet)
    {
        if(nDim<2)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cant' set Ny, dim=%3d !!!              ***\n",nDim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(ny<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ny=%6d is invalid                   ***\n",ny);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscFinalize();
        abort();
    }

    Ny=ny;
    IsNySet=true;
}
//**********************
void Mesh::SetNz(int nz)
{
    if(IsDimSet)
    {
        if(nDim<3)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cant' set Nz, dim=%3d !!!              ***\n",nDim);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
            PetscFinalize();
            abort();
        }
    }

    if(nz<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nz=%6d is invalid                   ***\n",nz);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
        PetscFinalize();
        abort();
    }

    Nz=nz;
    IsNzSet=true;
}

//*************************************
void Mesh::SetMeshType(string meshtype)
{
    if(IsDimSet)
    {
        if(nDim==1)
        {
            if(meshtype.compare(0,5,"edge2")==0)
            {
                MeshType="edge2";
                IsMeshTypeSet=true;
            }
            else if(meshtype.compare(0,5,"edge3")==0)
            {
                MeshType="edge3";
                IsMeshTypeSet=true;
            }
            else if(meshtype.compare(0,5,"edge4")==0)
            {
                MeshType="edge4";
                IsMeshTypeSet=true;
            }
            else
            {
                IsMeshTypeSet=false;
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported 1d mesh!                   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**         edge2,3,4 is valid for 1d mesh type!   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscFinalize();
                abort();
            }
        }
        else if(nDim==2)
        {
            if(meshtype.compare(0,5,"quad4")==0)
            {
                MeshType="quad4";
                IsMeshTypeSet=true;
            }
            else if(meshtype.compare(0,5,"quad8")==0)
            {
                MeshType="quad8";
                IsMeshTypeSet=true;
            }
            else if(meshtype.compare(0,5,"quad9")==0)
            {
                MeshType="quad9";
                IsMeshTypeSet=true;
            }
            else
            {
                IsMeshTypeSet=false;
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported 1d mesh!                   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**         quad4,8,9 is valid for 2d mesh type!   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscFinalize();
                abort();
            }
        }
        else if(nDim==3)
        {
            if(meshtype.compare(0,8,"hex8")==0)
            {
                MeshType="hex8";
                IsMeshTypeSet=true;
            }
            else if(meshtype.compare(0,5,"hex20")==0)
            {
                MeshType="hex20";
                IsMeshTypeSet=true;
            }
            else if(meshtype.compare(0,5,"hex27")==0)
            {
                MeshType="hex27";
                IsMeshTypeSet=true;
            }
            else
            {
                IsMeshTypeSet=false;
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported 1d mesh!                   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**         quad4,8,9 is valid for 2d mesh type!   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                   ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*****************************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else
    {
        MeshType=meshtype;
        IsMeshTypeSet=true;
    }
}