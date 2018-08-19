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
// get information of equation system in AsFem

#include "EquationSystem/EquationSystem.h"

int EquationSystem::GetDofIndexByName(string dofname) const
{
    bool FindName;
    if(!SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution has no name!!!         ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name to them   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    FindName=false;
    for(unsigned int i=0;i<solution_name_map.size();i++)
    {
        if(dofname==solution_name_map[i].first)
        {
            FindName=true;
            return solution_name_map[i].second;
        }
    }
    if(!FindName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: name=%10s can't be found!!        ***\n",dofname);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name to them   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return -1;
}


//*************************************************
string EquationSystem::GetDofNameByIndex(const int i) const
{
    bool FindIndex;
    if(!SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution has no name!!!         ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given name to them   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    FindIndex=false;
    for(unsigned int j=0;j<solution_name_map.size();j++)
    {
        if(solution_name_map[j].second==i)
        {
            FindIndex=true;
            return solution_name_map[j].first;
        }
    }

    if(!FindIndex)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: index=%3d can't be found!!      ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should set dofs index !!!   ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return "AsFem";
}
