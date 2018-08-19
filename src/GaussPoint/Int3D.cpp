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
// Created by walkandthinker on 19.08.18.
// Define standard 3D gauss point for 3d hex mesh

#include "GaussPoint/GaussPoint.h"

void Int3D(int ngp,int &Lint,double (&gs)[125][4])
{
    // generate gauss point for 3d integration
    if(ngp<2||ngp>6)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid gauss point number!!!   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Int3D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        abort();
    }

    //gs.setZero();
    Lint=ngp*ngp*ngp;


    int i,j,k,kk;

    double gs1d[125][4];

    Int1D(ngp,i,gs1d);
    kk=0;
    for(i=0;i<ngp;++i)
    {
        for(j=0;j<ngp;++j)
        {
            for(k=0;k<ngp;++k)
            {
                gs[kk][0]=gs1d[i][0]*gs1d[j][0]*gs1d[k][0];

                gs[kk][1]=gs1d[i][1];
                gs[kk][2]=gs1d[j][1];
                gs[kk][3]=gs1d[k][1];
                kk+=1;
            }
        }
    }
}


