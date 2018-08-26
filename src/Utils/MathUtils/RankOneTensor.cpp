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
// Created by walkandthinker on 26.08.18.
// Define the rank-1 tensor(vecotr) in AsFem

#include "Utils/RankOneTensor.h"

RankOneTensor::RankOneTensor(int ndim, double val)
{
    if(ndim<1||ndim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid in vector!   ***\n",ndim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    nDim=ndim;
    for(int i=0;i<nDim;i++) elements[i]=val;
}

//************************************
void RankOneTensor::SetDim(int dim)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid in vector!   ***\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    nDim=dim;
}

//*******************************
void RankOneTensor::ZeroEntities()
{
    for(int i=0;i<nDim;i++) elements[i]=0.0;
}

//******************************
double RankOneTensor::operator()(int i) const
{
    if(i<1||i>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range!          ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return elements[i-1];
}
//***********************************
double& RankOneTensor::operator()(int i)
{
    if(i<1||i>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range!          ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    return elements[i-1];
}
//*********************************
RankOneTensor RankOneTensor::operator+(double val)
{
    RankOneTensor temp(nDim,0.0);
    for(int i=0;i<nDim;i++)
        temp(i+1)=(*this)(i+1)+val;
    return temp;
}
RankOneTensor RankOneTensor::operator+(RankOneTensor &otherRankOneTensor)
{
    if(nDim!=otherRankOneTensor.GetSize())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied for two different vec!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    RankOneTensor temp(nDim,0.0);
    for(int i=0;i<nDim;i++)
        temp(i+1)=(*this)(i+1)+otherRankOneTensor(i+1);
    return temp;
}
//*********************************
RankOneTensor RankOneTensor::operator-(double val)
{
    RankOneTensor temp(nDim,0.0);
    for(int i=0;i<nDim;i++)
        temp(i+1)=(*this)(i+1)-val;
    return temp;
}
RankOneTensor RankOneTensor::operator-(RankOneTensor &otherRankOneTensor)
{
    if(nDim!=otherRankOneTensor.GetSize())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied for two different vec!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    RankOneTensor temp(nDim,0.0);
    for(int i=0;i<nDim;i++)
        temp(i+1)=(*this)(i+1)-otherRankOneTensor(i+1);
    return temp;
}
//*********************************
RankOneTensor RankOneTensor::operator*(double val)
{
    RankOneTensor temp(nDim,0.0);
    for(int i=0;i<nDim;i++)
        temp(i+1)=(*this)(i+1)*val;
    return temp;
}
double RankOneTensor::operator*(RankOneTensor &otherRankOneTensor)
{
    if(nDim!=otherRankOneTensor.GetSize())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: * applied for two different vec!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    double sum=0.0;
    for(int i=1;i<=nDim;i++)
    {
        sum+=(*this)(i)*otherRankOneTensor(i);
    }
    return sum;
}

//*******************************************
RankOneTensor& RankOneTensor::operator=(double val)
{
    for(int i=0;i<nDim;i++) elements[i]=val;
    return (*this);
}
RankOneTensor& RankOneTensor::operator=(RankOneTensor &otherRankOneTensor)
{
    if(nDim!=otherRankOneTensor.GetSize())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: = applied for two different vec!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=0;i<nDim;i++)
        elements[i]=otherRankOneTensor(i+1);

    return (*this);
}

//***************************************
void RankOneTensor::PrintTensor() const
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Rank-1 tensor:                         ***\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** ");
    for(int i=0;i<nDim;i++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%12.4e ",(*this)(i+1));
    }
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***\n");
}