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
// Created by walkandthinker on 18.08.18.
// uel for poisson equation

#include "ElementSystem/ElementSystem.h"

void ElementSystem::Poisson(const int &iState, const int (&IX)[27], const int &nDim, const int &nNodes,
                            const int &nDofs, const double &dt, const double &t, const double (&ctan)[2],
                            const double (&Coords)[27][4], const double (&U)[270][2], double (&K)[270][270],
                            double (&rhs)[270], double (&proj)[27][13])
{
    int i,j,k;
    if(iState%3==0)
    {
        // initializing local k and rhs
        for(i=0;i<nDofs;i++)
        {
            rhs[i]=0.0;
            for(j=0;j<nDofs;j++)
            {
                K[i][j]=0.0;
            }
        }
    }
}

