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
// define DofHanlder settings in AsFem

#include "DofHandler/DofHandler.h"

//*****************************
void DofHandler::SetDofsPerNode(int ndofspernode)
{
    if(ndofspernode<1||ndofspernode>10)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nDofsPerNode=%4d is invalid!!! ***\n",ndofspernode);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    nDofsPerNode=ndofspernode;

}

