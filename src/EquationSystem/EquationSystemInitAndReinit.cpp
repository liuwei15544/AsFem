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
// init and reinit funcions for equation system

#include "EquationSystem/EquationSystem.h"

PetscErrorCode EquationSystem::Init()
{
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    if(IsInit)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Warning: equation system is already initialized!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***          AsFem will do nothing for you!\n");

    }
    else
    {

        ierr=MatCreate(PETSC_COMM_WORLD,&AMATRIX);CHKERRQ(ierr);
        ierr=MatSetSizes(AMATRIX,PETSC_DECIDE,PETSC_DECIDE,nDofs,nDofs);CHKERRQ(ierr);
        NNZ_InExact=PetscInt(4*nDims*nNodesPerElmt*nDofsPerNode);// here I assume the each row
                                                 // has maximum 12 neighoured elmt,
                                                 // each elmt has maximum 27 node
        if(size>1)
        {
            MatSetType(AMATRIX,MATMPIAIJ);
            MatMPIAIJSetPreallocation(AMATRIX,NNZ_InExact,NULL,NNZ_InExact,NULL);
        }
        else
        {
            MatSetType(AMATRIX,MATSEQAIJ);
            MatSeqAIJSetPreallocation(AMATRIX,NNZ_InExact,NULL);
        }
        ierr=MatSetFromOptions(AMATRIX);CHKERRQ(ierr);
        ierr=MatSetUp(AMATRIX);CHKERRQ(ierr);

        // for vector

        VecCreate(PETSC_COMM_WORLD,&U0);//CHKERRQ(ierr);
        VecSetSizes(U0,PETSC_DECIDE,nDofs);//CHKERRQ(ierr);
        VecSetFromOptions(U0);//CHKERRQ(ierr);
        VecSetUp(U0);//CHKERRQ(ierr);


        VecDuplicate(U0,&RHS);//CHKERRQ(ierr);
        VecDuplicate(U0,&dU);//CHKERRQ(ierr);
        VecDuplicate(U0,&U);//CHKERRQ(ierr);
        VecDuplicate(U0,&V);//CHKERRQ(ierr);

        // initialize
        const PetscScalar val=1.0e-12;
        ierr=VecSet(U0,val);//CHKERRQ(ierr);
        ierr=VecSet(RHS,val);//CHKERRQ(ierr);
        ierr=VecSet(V,val);//CHKERRQ(ierr);



        nDofsPerNode=int(solution_name_map.size());
        int n=nDofs/nDofsPerNode;

        VecCreate(PETSC_COMM_WORLD,&Proj);//CHKERRQ(ierr);
        VecSetSizes(Proj,PETSC_DECIDE,n*13);//CHKERRQ(ierr);
        VecSetFromOptions(Proj);//CHKERRQ(ierr);
        VecSetUp(Proj);//CHKERRQ(ierr);
        VecSet(Proj,val);

        return 1;
    }

    return ierr;
}

//**********************************************
void EquationSystem::ReInitKandR()
{
    // Warning: here this function should and must only
    //          initializing the AMATRIX and RHS,
    //          don't do operating on other vec!!!

    // so this should be called before the N-R iteration!

    MatZeroEntries(AMATRIX);// must be sure AMATRIX is already initalized
    VecSet(RHS,0.0);
}

//*******************************************
void EquationSystem::ReInitVec(Vec &v)
{
    VecSet(v,0.0);
}
//****************************************
void EquationSystem::ReInitMat(Mat &a)
{
    MatZeroEntries(a);// must be sure AMATRIX is already initalized
}

//****************************************
void EquationSystem::UpdateUplusdU()
{
    VecAXPY(U,1.0,dU);
}
//****************************************
void EquationSystem::UpdateU0(Vec &u)
{
    //VecCopy(Vec x, Vec y); y = x
    VecCopy(u,U0);
}
//*****************************************
void EquationSystem::UpdateV(Vec &v)
{
    VecCopy(v,V);
}



