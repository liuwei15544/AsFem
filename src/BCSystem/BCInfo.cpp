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
// define the boundary info for input file(combine all the bc block)

#include "BCSystem/BCInfo.h"

BCInfo::BCInfo()
{
    HasBCKernelBlocks=false;
    BCBlockDofIndex.clear();
    BCBlockList.clear();
}

//*************************
void BCInfo::Init()
{
    HasBCKernelBlocks=false;
    BCBlockDofIndex.clear();
    BCBlockList.clear();
}

//*************************
void BCInfo::AddBCKernelBlock(BCBlockInfo &bcBlockInfo)
{
    if(BCBlockList.size()==0)
    {
        BCBlockList.push_back(bcBlockInfo);
        HasBCKernelBlocks=true;
    }
    else
    {
        bool IsBlockNameInList=false;
        for(size_t i=0;i<BCBlockList.size();i++)
        {
            if(bcBlockInfo.BCBlockName==BCBlockList[i].BCBlockName)
            {
                IsBlockNameInList=true;
                break;
            }
        }
        if(!IsBlockNameInList)
        {
            BCBlockList.push_back(bcBlockInfo);
            HasBCKernelBlocks=true;
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: duplicated bc block info!!!     ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
}

void BCInfo::AddBCKernelBlockFromList(vector<BCBlockInfo> bcBlockList)
{
    for(unsigned int i=0;i<bcBlockList.size();i++)
    {
        AddBCKernelBlock(bcBlockList[i]);
    }
}

//*********************************
BCBlockInfo BCInfo::GetIthBCKernelBlock(int i) const
{
    if(!HasBCKernelBlocks)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: no NodalBCKernelBlocks          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>int(BCBlockList.size()))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range!!!        ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    else
    {
        return BCBlockList[i-1];
    }
}


bool BCInfo::CheckDofsName(EquationSystem &equationSystem)
{
    DofsNameCheckPassed=false;
    if(!HasBCKernelBlocks)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: BCKernelBlockList is empty !!!  ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't check dofsname of bckernel***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    bool IsDofsNameValid=false;
    string dof;
    for(unsigned int i=0;i<BCBlockList.size();i++)
    {
        dof=BCBlockList[i].BCBlockDofsName;
        IsDofsNameValid=false;
        for(int k=0;k<equationSystem.GetDofsNumPerNode();k++)
        {
            if(dof==equationSystem.GetDofNameByIndex(k+1))
            {
                IsDofsNameValid=true;
                break;
            }
        }
        if(!IsDofsNameValid)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: in %3d-th bc kernel block!!     ***\n",i+1);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dof name is not correct!!!      ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }

    return IsDofsNameValid;

}
//**************************
void BCInfo::GenerateBCKernelDofMap(EquationSystem &equationSystem)
{
    if(CheckDofsName(equationSystem))
    {
        BCBlockDofIndex.clear();
        int i,j,k,iInd;
        string name;
        bool IsNameInList=false;
        for(i=0;i<BCBlockList.size();i++)
        {
            name=BCBlockList[i].BCBlockDofsName;
            IsNameInList=false;
            for(j=0;j<equationSystem.GetDofsNumPerNode();j++)
            {
                if(name==equationSystem.GetDofNameByIndex(j+1))
                {
                    iInd=equationSystem.GetDofIndexByName(name);
                    IsNameInList=true;
                    break;
                }
            }
            if(!IsNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: in %3d-th bc kernel block!!     ***\n",i+1);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dof name is not correct!!!      ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                BCBlockDofIndex.push_back(iInd);
            }
        }
    }
}
//*************************
void BCInfo::PrintBCInfo() const
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Boundary block information:            ***\n");
    for(size_t i=0;i<BCBlockList.size();i++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   block name =%16s         ***\n",BCBlockList[i].BCBlockName.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***     kernel name=%16s       ***\n",BCBlockList[i].BCBlockKernelName.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***     side name  =%16s       ***\n",BCBlockList[i].sidename.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***     applied dof=%16s       ***\n",BCBlockList[i].BCBlockDofsName.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
}

//***********************************************
int BCInfo::GetIthBCKernelDofIndex(int i) const
{
    if(i<1||i>BCBlockList.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of bc kernel range!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return BCBlockDofIndex[i-1];
}

//***************************
string BCInfo::GetIthBCKernelName(int i) const
{
    if(i<1||i>BCBlockList.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of bc kernel range!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return BCBlockList[i-1].BCBlockKernelName;
}

//*************************
double BCInfo::GetIthBCKernelValue(int i) const
{
    if(i<1||i>BCBlockList.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of bc kernel range!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return BCBlockList[i-1].BCValue;
}

//***********************
string BCInfo::GetIthBCKernelSideName(int i) const
{
    if(i<1||i>BCBlockList.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of bc kernel range!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return BCBlockList[i-1].sidename;
}


