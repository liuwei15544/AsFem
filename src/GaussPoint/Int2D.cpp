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
// Define standard 2D gauss point for 2d quadratic mesh

#include "GaussPoint/GaussPoint.h"

void Int2D(int ngp,int &Lint,double (&gs)[125][4])
{
    // generate gauss point for 2d integration
    if(ngp<2||ngp>6)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid gauss point number!!!   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Int2D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        abort();
    }

    Lint=ngp*ngp;

    if (ngp == 2)
    {
        // 3----4
        // |    |
        // 1----2
        gs[0][1]=-0.57735026918962573;
        gs[0][2]=-0.57735026918962573;
        gs[0][0]= 1.0;

        gs[1][1]= 0.57735026918962573;
        gs[1][2]=-0.57735026918962573;
        gs[1][0]= 1.0;

        gs[2][1]= 0.57735026918962573;
        gs[2][2]= 0.57735026918962573;
        gs[2][0]= 1.0;

        gs[3][1]=-0.57735026918962573;
        gs[3][2]= 0.57735026918962573;
        gs[3][0]= 1.0;
    }
    else if (ngp==3)
    {
        // 7--8--9
        // |  |  |
        // 4--5--6
        // |  |  |
        // 1--2--3
        double w1=5.0/9.0,w2=8.0/9.0;

        gs[0][1]=-0.7745966692414834;
        gs[0][2]=-0.7745966692414834;
        gs[0][0]= w1*w1;

        gs[1][1]= 0.7745966692414834;
        gs[1][2]=-0.7745966692414834;
        gs[1][0]= w1*w1;

        gs[2][1]= 0.7745966692414834;
        gs[2][2]= 0.7745966692414834;
        gs[2][0]= w1*w1;

        gs[3][1]=-0.7745966692414834;
        gs[3][2]= 0.7745966692414834;
        gs[3][0]= w1*w1;

        gs[4][1]= 0.0;
        gs[4][2]=-0.7745966692414834;
        gs[4][0]= w1*w2;

        gs[5][1]= 0.7745966692414834;
        gs[5][2]= 0.0;
        gs[5][0]= w1*w2;

        gs[6][1] = 0.0;
        gs[6][2] = 0.7745966692414834;
        gs[6][0] = w1*w2;

        gs[7][1] = -0.7745966692414834;
        gs[7][2] = 0.0;
        gs[7][0] = w1*w2;

        gs[8][1] = 0.0;
        gs[8][2] = 0.0;
        gs[8][0] = w2*w2;
    }
    else if(ngp>=4)
    {
        double gs1d[125][4];

        int i,j,k;
        Int1D(ngp,i,gs1d);
        k=0;
        for(j=0;j<ngp;++j)
        {
            for(i=0;i<ngp;++i)
            {
                gs[k][0]=gs1d[i][0]*gs1d[j][0];
                gs[k][1]=gs1d[i][1];
                gs[k][2]=gs1d[j][1];
                k+=1;
            }
        }
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid gaus point number!!!    ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Int1D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        abort();
    }
}

