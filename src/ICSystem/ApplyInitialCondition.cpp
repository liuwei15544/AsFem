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
// Created by walkandthinker on 25.08.18.
// apply initial condition system in AsFem

#include "ICSystem/ICSystem.h"

void ICSystem::ApplyInitialCondition(Mesh &mesh, Vec &U)
{
    int i;
    for(i=0;i<icInfo.GetICBlockNum();i++)
    {
        SingleICBlock=icInfo.GetIthICKernelBlock(i+1);
        if(icInfo.GetIthICKernelName(i+1)=="constic")
        {
            ConstantIC(SingleICBlock,mesh,U);
        }
        else if(icInfo.GetIthICKernelName(i+1)=="circleic")
        {
            CircleIC(SingleICBlock,mesh,U);
        }
        else if(icInfo.GetIthICKernelName(i+1)=="randomic")
        {
            RandomIC(SingleICBlock,mesh,U);
        }
        else if(icInfo.GetIthICKernelName(i+1)=="randomnoiseic")
        {
            RandomNoiseIC(SingleICBlock,mesh,U);
        }
    }


    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
}

