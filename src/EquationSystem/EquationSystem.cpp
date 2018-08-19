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
// Define equation system in AsFem

#include "EquationSystem/EquationSystem.h"

EquationSystem::EquationSystem()
{

    nDofs=0;nDofsPerNode=0;
    solution_name_list.clear();
    solution_index_list.clear();
    IsInit=false;
    SolutionHasName=false;
}

//***************************************
EquationSystem::EquationSystem(const int dofs,int dofspernode)
{

    if(dofs<2||dofs>MaxDofs)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dofs number=%-8d is not supported in current version!\n",dofs);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        the maximum dofs of current version is:%-8d!\n",MaxDofs);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(dofspernode<1||dofspernode>10)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nDofsPerNode=%d is invalid in current version!!!\n",dofspernode);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        max dofs per node<=10!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    nDofs=dofs;nDofsPerNode=dofspernode;
    solution_name_list.clear();
    solution_index_list.clear();
    IsInit=false;
    SolutionHasName=false;
}


//************************************************
PetscErrorCode EquationSystem::Release()
{
    ierr=MatDestroy(&AMATRIX);CHKERRQ(ierr);

    ierr=VecDestroy(&U0);CHKERRQ(ierr);
    ierr=VecDestroy(&U);CHKERRQ(ierr);
    ierr=VecDestroy(&V);CHKERRQ(ierr);
    ierr=VecDestroy(&dU);CHKERRQ(ierr);
    ierr=VecDestroy(&RHS);CHKERRQ(ierr);

    solution_name_list.clear();
    solution_index_list.clear();

    return ierr;
}

//**************************************************
//*** Print out the info of equation system
//**************************************************
void EquationSystem::PrintSolutionNameMap(string str) const
{
    if(!SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution has no name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should set name for it first!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    if(str.size()>1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** str= %s\n",str.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Solution system information:           ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   number of solution=%2d                ***\n",nDofsPerNode);

    string name;
    for(int i=1;i<=nDofsPerNode;i++)
    {
        name=GetDofNameByIndex(i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   solution[%2d]<------>%8s         ***\n",i,name.c_str());
    }

    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");

}