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
// define the projection function from gauss point to nodal point

#include "ElementSystem/ElementSystem.h"

void ElementSystem::Projection(const int &nNodes, const double &xsj, double (&shp)[27][4], const double (&value)[12],
                               double (&proj)[27][13])
{
    int i,k;
    double weight;
    for(i=0;i<nNodes;i++)
    {
        weight=shp[i][0]*xsj;
        proj[i][0]+=weight;
        for(k=1;k<=12;k++)
        {
            proj[i][k]+=value[k-1]*weight;
        }
    }
}

