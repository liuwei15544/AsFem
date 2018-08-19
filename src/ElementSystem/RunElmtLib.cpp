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
// the element library in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::RunElmtLib(const int &iState, const int (&IX)[27], const int &nDim, const int &nNodes,
                               const int &nDofs, const double &dt, const double &t, const double (&ctan)[2],
                               const double (&Coords)[27][4], const double (&U)[270][2], double (&K)[270][270],
                               double (&rhs)[270], double (&proj)[27][12+1])
{
    switch (ActiveUelIndex)
    {
        case UelList::poisson:
            Poisson(iState,IX,nDim,nNodes,nDofs,dt,t,ctan,Coords,U,K,rhs,proj);
            break;
        case UelList::diffusion:
            break;
        case UelList::cahnhilliard:
            break;
        case UelList::mechanics:
            break;
        case UelList::thermalmechanics:
            break;
        case UelList::uel1:
            break;
        case UelList::uel2:
            break;
        case UelList::uel3:
            break;
        case UelList::uel4:
            break;
        case UelList::uel5:
            break;
        case UelList::uel6:
            break;
        case UelList::uel7:
            break;
        case UelList::uel8:
            break;
        case UelList::uel9:
            break;
        case UelList::uel10:
            break;
        default:
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported uel in ElmtLib!!!   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
    }
}

