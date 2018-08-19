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
// print out dof information in AsFem

#include "DofHandler/DofHandler.h"

void DofHandler::PrintDofMap() const
{
    if(!HasDofMap)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't print dof map info, DofHandler hasn't been generated!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should create mesh first, then generate the dof map\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        then AsFem can print dof map for you!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Dof map information:                   ***\n");
    for(int e=1;e<=nElmts;e++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** e=%6d: ",e);
        for(int i=1;i<=nDofsPerElmt;i++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",GlobalDofMap[(e-1)*nDofsPerElmt+i-1]);
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Boundary dof map information:          ***\n");

    int dofind[270],len;
    string sidename;
    for(unsigned int i=0;i<GlobalBCDofMap.size();i++)
    {
        sidename=GlobalBCDofMap[i].first;
        for(int e=1;e<=GetBCSideDofsNum(sidename)/nDofsPerBCElmt;e++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** %8s side->",sidename.c_str());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"e=%3d:",e);

            GetLocalBCDofMap(sidename,e,len,dofind);
            for(int j=0;j<len;j++)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%5d ",dofind[j]);
            }
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
        }
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
}



