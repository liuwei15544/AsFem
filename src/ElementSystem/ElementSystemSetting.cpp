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
// define settings for the uel system in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::SetUelIndex()
{
    if(kernelBlockInfo.ElementName=="poisson")
    {
        ActiveUelIndex=UelList::poisson;
    }
    else if(kernelBlockInfo.ElementName=="diffusion")
    {
        ActiveUelIndex=UelList::diffusion;
    }
    else if(kernelBlockInfo.ElementName=="cahnhilliard")
    {
        ActiveUelIndex=UelList::cahnhilliard;
    }
    else if(kernelBlockInfo.ElementName=="allencahn")
    {
        ActiveUelIndex=UelList::allencahn;
    }
    else if(kernelBlockInfo.ElementName=="mechanics")
    {
        ActiveUelIndex=UelList::mechanics;
    }
    else if(kernelBlockInfo.ElementName=="thermalmechanics")
    {
        ActiveUelIndex=UelList::thermalmechanics;
    }
    else if(kernelBlockInfo.ElementName=="uel1")
    {
        ActiveUelIndex=UelList::uel1;
    }
    else if(kernelBlockInfo.ElementName=="uel2")
    {
        ActiveUelIndex=UelList::uel2;
    }
    else if(kernelBlockInfo.ElementName=="uel3")
    {
        ActiveUelIndex=UelList::uel3;
    }
    else if(kernelBlockInfo.ElementName=="uel4")
    {
        ActiveUelIndex=UelList::uel4;
    }
    else if(kernelBlockInfo.ElementName=="uel5")
    {
        ActiveUelIndex=UelList::uel5;
    }
    else if(kernelBlockInfo.ElementName=="uel6")
    {
        ActiveUelIndex=UelList::uel6;
    }
    else if(kernelBlockInfo.ElementName=="uel7")
    {
        ActiveUelIndex=UelList::uel7;
    }
    else if(kernelBlockInfo.ElementName=="uel8")
    {
        ActiveUelIndex=UelList::uel8;
    }else if(kernelBlockInfo.ElementName=="uel9")
    {
        ActiveUelIndex=UelList::uel9;
    }else if(kernelBlockInfo.ElementName=="uel10")
    {
        ActiveUelIndex=UelList::uel10;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported uel type!!!         ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
}

//************************************************************
void ElementSystem::SetUmatIndex()
{
    if(kernelBlockInfo.MaterialKernelName.size()<1)
    {
        ActiveUmatIndex=UmatList::emptyumat;
    }
    else if(kernelBlockInfo.MaterialKernelName=="linearelastic")
    {
        ActiveUmatIndex=UmatList::linearelastic;
    }
    else if(kernelBlockInfo.MaterialKernelName=="deformgrad")
    {
        ActiveUmatIndex=UmatList::deformgrad;
    }
    else if(kernelBlockInfo.MaterialKernelName=="freeenergy")
    {
        ActiveUmatIndex=UmatList::freeenergy;
    }
    else if(kernelBlockInfo.MaterialKernelName=="freeenergyac")
    {
        ActiveUmatIndex=UmatList::freeenergyac;
    }
    else if(kernelBlockInfo.MaterialKernelName=="conductivity")
    {
        ActiveUmatIndex=UmatList::conductivity;
    }
    else if(kernelBlockInfo.MaterialKernelName=="neohookean")
    {
        ActiveUmatIndex=UmatList::neohookean;
    }
    else
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported umat type!!!        ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(kernelBlockInfo.strain=="small")
    {
        strainMode=small;
    }
    else
    {
        strainMode=finite;
    }
}