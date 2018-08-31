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
// Define 1d shape function used in AsFem

#include "ShapeFuns/ShapeFuns.h"

void Shp1D(const int &ndim,const int &nnodes,const double &xi,const double (&Coords)[27][4],
           double (&shp)[27][4],double &DetJac)
{
    double dxdxi,dydxi,dzdxi;
    int i;

    if (ndim<1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Wrong dimension case !!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Shp1D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");

        abort();
    }

    if (nnodes<2 || nnodes>4)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Wrong dimension case !!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Wrong nodes num per element(nNodes=%2d)!!!\n",nnodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Shp1D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem only support 1-3 order 1D mesh!!!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        abort();
    }

    //shp.setZero();
    if (nnodes == 2)
    {
        // Linear element
        shp[0][0] = 0.5*(1.0 - xi);
        shp[0][1] =-0.5;

        shp[1][0] = 0.5*(xi + 1.0);
        shp[1][1] = 0.5;
    }
    else if (nnodes == 3)
    {
        // Quadratic line element(3 nodes)
        shp[0][0] = 0.5*xi*(xi - 1.0);
        shp[0][1] = 0.5*(2.0*xi - 1.0);

        shp[1][0] = -(xi + 1.0)*(xi - 1.0);
        shp[1][1] = -2.0*xi;

        shp[2][0] = 0.5*xi*(xi + 1.0);
        shp[2][1] = 0.5*(2.0*xi + 1.0);
    }
    else if(nnodes==4)
    {
        // Third order line element
        shp[0][0]=-(3.0*xi+1.0)*(3.0*xi-1.0)*(    xi-1.0)/16.0;
        shp[1][0]= (3.0*xi+3.0)*(3.0*xi-1.0)*(3.0*xi-3.0)/16.0;
        shp[2][0]=-(3.0*xi+3.0)*(3.0*xi+1.0)*(3.0*xi-3.0)/16.0;
        shp[3][0]= (    xi+1.0)*(3.0*xi+1.0)*(3.0*xi-1.0)/16.0;
        shp[0][1]=-27.0*xi*xi/16.0+9.0*xi/8.0+ 1.0/16.0;
        shp[1][1]= 81.0*xi*xi/16.0-9.0*xi/8.0-27.0/16.0;
        shp[2][1]=-81.0*xi*xi/16.0-9.0*xi/8.0+27.0/16.0;
        shp[3][1]= 27.0*xi*xi/16.0+9.0*xi/8.0- 1.0/16.0;
    }


    // Compute the derivatives over x
    // dN/dx
    // actually, there should be positive and negative
    // but I want user to decide the direction of normal vector
    // so the DetJac is always positive here
    dxdxi=0.0;dydxi=0.0;dzdxi=0.0;
    for (i=0;i<nnodes;++i)
    {
        dxdxi+=Coords[i][1]*shp[i][1];
        dydxi+=Coords[i][2]*shp[i][1];
        dzdxi+=Coords[i][3]*shp[i][1];
    }
    DetJac=sqrt(dxdxi*dxdxi+dydxi*dydxi+dzdxi*dzdxi);


    if (fabs(DetJac)<1.e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: Singular in element(Shp1D)!!!   ***\n");
        abort();
    }


    for (i=0;i<nnodes;++i)
    {
        shp[i][1]=shp[i][1]/DetJac;
    }
}

