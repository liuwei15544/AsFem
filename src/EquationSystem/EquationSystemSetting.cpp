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
// Define settings for equation system in AsFem

#include "EquationSystem/EquationSystem.h"

void EquationSystem::AddSolutionNameAndIndex(string name, int order)
{
    if(SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution already has name!!!\n",name,order);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't add name and order to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    bool IsNameIn=false,IsOrderIn=false;
    for(unsigned int i=0;i<solution_name_list.size();i++)
    {
        if(name==solution_name_list[i])
        {
            IsNameIn=true;
            break;
        }
    }

    for(unsigned int i=0;i<solution_index_list.size();i++)
    {
        if(order==solution_index_list[i])
        {
            IsOrderIn=true;
            break;
        }
    }

    if(IsNameIn || IsOrderIn)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: name=%8s or order=%2d is already in the list!!\n",name.c_str(),order);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        unique name and oder should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    else
    {
        solution_name_list.push_back(name);
        solution_index_list.push_back(order);
    }
}

void EquationSystem::SetSolutionNameFromVector(vector<string> names, vector<int> orders)
{
    if(SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution already has name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't add name and order to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(names.size()!=orders.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: name list isnt't match oder list\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        same size vector should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(nDofsPerNode==0) nDofsPerNode=int(names.size());


    bool IsNameUnique,IsOrderUnique;

    auto its=unique(names.begin(),names.end());
    IsNameUnique=(its==names.end());

    auto it=unique(orders.begin(),orders.end());
    IsOrderUnique=(it==orders.end());

    if((!IsNameUnique) || (!IsOrderUnique))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: either name or oder list isn't unique\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        unique name and oder should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    if(names.size()!=orders.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: len(names)=%2d isn't equal len(orders)=%2d\n",names.size(),orders.size());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        equal length names_list and orders_list should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    int ii=int(names.size());
    int jj=int(orders.size());
    if(ii<nDofsPerNode || jj<nDofsPerNode)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: neighter name or oder list is long enough to set solution name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please add more name and order values to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }
    else
    {

        pair<string,int> temp;
        for(unsigned int i=0;i<orders.size();i++)
        {
            temp=make_pair(names[i],orders[i]);
            solution_name_map.push_back(temp);
        }
        SolutionHasName=true;
        nDofsPerNode=orders.size();
    }
}

void EquationSystem::SetSolutionName()
{
    if(SolutionHasName)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: solution already has name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't add name and order to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    int ii=int(solution_index_list.size());
    int jj=int(solution_name_list.size());
    if(ii<nDofsPerNode || jj<nDofsPerNode)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: neighter name or oder list is long enough to set solution name!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please add more name and order values to the list!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    bool IsNameUnique,IsOrderUnique;

    auto its=unique(solution_name_list.begin(),solution_name_list.end());
    IsNameUnique=(its==solution_name_list.end());

    auto it=unique(solution_index_list.begin(),solution_index_list.end());
    IsOrderUnique=(it==solution_index_list.end());

    if((!IsNameUnique) || (!IsOrderUnique))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: either name or oder list isn't unique\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        unique name and oder should be given!!!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!\n");
        PetscFinalize();
        abort();
    }

    pair<string,int> temp;
    for(unsigned int i=0;i<solution_index_list.size();i++)
    {
        temp=make_pair(solution_name_list[i],solution_index_list[i]);
        solution_name_map.push_back(temp);
    }
    SolutionHasName=true;
    nDofsPerNode=int(solution_index_list.size());
}


