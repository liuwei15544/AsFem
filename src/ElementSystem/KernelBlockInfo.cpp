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
// define the kernel block for input file reader

#include "ElementSystem/KernelBlockInfo.h"

KernelBlockInfo::KernelBlockInfo()
{
    ElementName="";
    MaterialKernelName="";
    params.clear();
}
//*********************************
void KernelBlockInfo::Init()
{
    ElementName="";
    MaterialKernelName="";
    params.clear();
}
//*********************************
void KernelBlockInfo::PrintKernelBlockInfo() const
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Kernel block information:              ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   elmt  name= %-18s       ***\n",ElementName.c_str());
    if(MaterialKernelName.size()>0)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   material name= %-18s    ***\n",MaterialKernelName.c_str());
    }
    if(params.size()>0)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   params=");
        for(unsigned int i=0;i<params.size();i++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%8.3e ",params[i]);
            if(params.size()>5)
            {
                if((i+1)%5==0)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n             ");
                }
            }
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");

    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
}


