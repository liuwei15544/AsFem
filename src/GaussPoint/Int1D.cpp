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
// Define standard 1D gauss point for 1d lagrange mesh

#include "GaussPoint/GaussPoint.h"

void Int1D(int ngp,int &Lint,double (&gs)[125][4])
{
    // generate gauss point for 1d integration
    if (ngp < 2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid gaus point number!!!    ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Int1D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        abort();
    }

    //gs.setZero();
    Lint=ngp;

    if (ngp == 2)
    {
        gs[0][1] =-0.577350269189625764509148780502;
        gs[1][1] = 0.577350269189625764509148780502;

        gs[0][0] = 1.0;
        gs[1][0] = 1.0;
    }
    else if (ngp == 3)
    {
        gs[0][0] = 0.555555555555555555555555555556;
        gs[0][1] =-0.774596669241483377035853079956;

        gs[1][0] = 0.888888888888888888888888888889;
        gs[1][1] = 0.0;

        gs[2][0] = 0.555555555555555555555555555556;
        gs[2][1] = 0.774596669241483377035853079956;

    }
    else if (ngp == 4)
    {
        gs[0][0] = 0.347854845137453857373063949222;
        gs[0][1] =-0.861136311594052575223946488893;

        gs[1][0]= 0.652145154862546142626936050778;
        gs[1][1]=-0.339981043584856264802665759103;

        gs[2][0]= 0.652145154862546142626936050778;
        gs[2][1]= 0.339981043584856264802665759103;

        gs[3][0]= 0.347854845137453857373063949222;
        gs[3][1]= 0.861136311594052575223946488893;
    }
    else if (ngp==5)
    {
            gs[0][1]=-0.906179845938663992797626878299;
            gs[1][1]=-0.538469310105683091036314420700;
            gs[2][1]= 0.000000000000000000000000000000;
            gs[3][1]= 0.538469310105683091036314420700;
            gs[4][1]= 0.906179845938663992797626878299;

            gs[0][0]= 0.236926885056189087514264040720;
            gs[1][0]= 0.478628670499366468041291514836;
            gs[2][0]= 0.568888888888888888888888888889;
            gs[3][0]= 0.478628670499366468041291514836;
            gs[4][0]= 0.236926885056189087514264040720;
    }
//        else if(ngp==6)
//        {
//            gs.coeffRef(1,0)=-0.932469514203152027812301554494;
//            gs.coeffRef(1,1)=-0.661209386466264513661399595020;
//            gs.coeffRef(1,2)=-0.238619186083196908630501721681;
//            gs.coeffRef(1,3)= 0.238619186083196908630501721681;
//            gs.coeffRef(1,4)= 0.661209386466264513661399595020;
//            gs.coeffRef(1,5)= 0.932469514203152027812301554494;
//
//            gs.coeffRef(0,0)= 0.171324492379170345040296142173;
//            gs.coeffRef(0,1)= 0.360761573048138607569833513838;
//            gs.coeffRef(0,2)= 0.467913934572691047389870343990;
//            gs.coeffRef(0,3)= 0.467913934572691047389870343990;
//            gs.coeffRef(0,4)= 0.360761573048138607569833513838;
//            gs.coeffRef(0,5)= 0.171324492379170345040296142173;
//        }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: npg>5 is not supported !!!      ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        error happen in Int1D!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        abort();
    }
}

