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
// Created by walkandthinker on 25.08.18.
// define the initial condition info for AsFem

#include "ICSystem/ICInfo.h"

ICInfo::ICInfo()
{
    HasICKernelBlocks=false;
    ICBlockDofIndex.clear();
    ICBlockList.clear();
}
//****************
void ICInfo::Init()
{
    HasICKernelBlocks=false;
    ICBlockDofIndex.clear();
    ICBlockList.clear();
}

//**************************
void ICInfo::AddICKernelBlock(ICBlockInfo &icBlockInfo)
{
    if(ICBlockList.size()==0)
    {
        ICBlockList.push_back(icBlockInfo);
        HasICKernelBlocks=true;
    }
    else
    {
        bool IsBlockNameInList=false;
        for(unsigned int i=0;i<ICBlockList.size();i++)
        {
            if(icBlockInfo.ICBlockName==ICBlockList[i].ICBlockName)
            {
                IsBlockNameInList=true;
                break;
            }
        }
        if(!IsBlockNameInList)
        {
            ICBlockList.push_back(icBlockInfo);
            HasICKernelBlocks=true;
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: duplicated ic block info!!!     ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
}

//**************************************************
void ICInfo::AddICKernelBlockFromList(vector<ICBlockInfo> icBlockList)
{
    for(unsigned int i=0;i<icBlockList.size();i++)
    {
        AddICKernelBlock(icBlockList[i]);
    }
}
//**************************************************
ICBlockInfo ICInfo::GetIthICKernelBlock(int i) const
{
    if(!HasICKernelBlocks)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: no ICKernelBlock found!         ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>int(ICBlockList.size()))
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
        return ICBlockList[i-1];
    }
}
//****************************
bool ICInfo::CheckDofsName(EquationSystem &equationSystem)
{
    DofsNameCheckPassed=false;
    if(!HasICKernelBlocks)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ICKernelBlockList is empty !!!  ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't check dofsname of ickernel***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    bool IsDofsNameValid=false;
    string dof;
    for(unsigned int i=0;i<ICBlockList.size();i++)
    {
        dof=ICBlockList[i].ICBlockDofsName;
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
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: in %3d-th ic kernel block!!     ***\n",i+1);
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
void ICInfo::GenerateICKernelDofMap(EquationSystem &equationSystem)
{
    if(CheckDofsName(equationSystem))
    {
        ICBlockDofIndex.clear();
        int i,j,iInd;
        string name;
        bool IsNameInList=false;
        for(i=0;i<int(ICBlockList.size());i++)
        {
            name=ICBlockList[i].ICBlockDofsName;
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
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: in %3d-th ic kernel block!!     ***\n",i+1);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        dof name is not correct!!!      ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                ICBlockDofIndex.push_back(iInd);
                ICBlockList[i].ICBlockDofsIndex=iInd;
            }
        }
    }
}
//*************************
void ICInfo::PrintICInfo() const
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Boundary block information:            ***\n");
    for(unsigned int i=0;i<ICBlockList.size();i++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   block name =%16s         ***\n",ICBlockList[i].ICBlockName.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***     kernel name=%16s       ***\n",ICBlockList[i].ICBlockKernelName.c_str());
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***     applied dof=%16s       ***\n",ICBlockList[i].ICBlockDofsName.c_str());
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
}

//***********************************************
int ICInfo::GetIthICKernelDofIndex(int i) const
{

    if(i<1||i>int(ICBlockList.size()))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of ic kernel range!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    //TODO: Error , size of BCBlockDofIndex=0

    return ICBlockDofIndex[i-1];
}

//***************************
string ICInfo::GetIthICKernelName(int i) const
{
    if(i<1||i>int(ICBlockList.size()))
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of ic kernel range!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return ICBlockList[i-1].ICBlockKernelName;
}

//*************************
//double ICInfo::GetIthICKernelValue(int i) const
//{
//    if(i<1||i>int(ICBlockList.size()))
//    {
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of ic kernel range!***\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
//        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
//        PetscFinalize();
//        abort();
//    }
//    return BCBlockList[i-1].BCValue;
//}


