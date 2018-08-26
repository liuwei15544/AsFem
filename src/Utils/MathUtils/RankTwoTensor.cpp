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
// Define the rank-2 tensor(i.e. stress and strain) in AsFem

#include "Utils/RankTwoTensor.h"

//****************************************************
//*** For constructor
//****************************************************
RankTwoTensor::RankTwoTensor(int dim, double val)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid in rank2 tensor*\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    nDim=dim;
    for(int i=0;i<9;i++) elements[i]=val;
}
//***
RankTwoTensor::RankTwoTensor(int dim,RankTwoTensor::InitMethod method)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid in rank2 tensor*\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    nDim=dim;
    if(method==InitZero)
    {
        for(int i=0;i<9;i++) elements[i]=0.0;
    }
    else if(method==InitIdentity)
    {
        for(int i=0;i<nDim;i++)
        {
            for(int j=0;j<nDim;j++)
            {
                if(i==j)
                {
                    elements[i*nDim+j]=1.0;
                }
                else
                {
                    elements[i*nDim+j]=0.0;
                }
            }
        }
    }
}
//*** Construct from voigt notation
RankTwoTensor::RankTwoTensor(double &s11, double &s22, double &s12)
{
    // For 2D case
    nDim=2;
    (*this)(1,1)=s11;(*this)(1,2)=s12;
    (*this)(2,1)=s12;(*this)(2,2)=s22;

}
RankTwoTensor::RankTwoTensor(double &s11,double &s22,double &s33,double &s23,double &s13,double &s12)
{
    // For 3D case
    nDim=3;
    (*this)(1,1)=s11;(*this)(1,2)=s12;(*this)(1,3)=s13;
    (*this)(2,1)=s12;(*this)(2,2)=s22;(*this)(2,3)=s23;
    (*this)(3,1)=s13;(*this)(3,2)=s23;(*this)(3,3)=s33;
}
//*** Construct from displacement gradient(or other gradient, here mainly for the F tensor)
RankTwoTensor::RankTwoTensor(double (&gradUx)[3],double (&gradUy)[3])
{
    // For 2D displacement gradient
    nDim=2;
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];
}
RankTwoTensor::RankTwoTensor(double (&gradUx)[3],double (&gradUy)[3],double (&gradUz)[3])
{
    // For 3D displacement gradient
    nDim=3;
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];(*this)(1,3)=gradUx[2];
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];(*this)(2,3)=gradUy[2];
    (*this)(3,1)=gradUz[0];(*this)(3,2)=gradUz[1];(*this)(3,3)=gradUz[2];
}
//******************************************************
//*** Operator overload
//******************************************************
double RankTwoTensor::operator()(int i, int j) const
{
    if(i<1||i>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range in rank2  ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(j<1||j>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%3d is out of range in rank2  ***\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return elements[(i-1)*nDim+j-1];
}
double& RankTwoTensor::operator()(int i, int j)
{
    if(i<1||i>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range in rank2  ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(j<1||j>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%3d is out of range in rank2  ***\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return elements[(i-1)*nDim+j-1];
}
//********************************************
RankTwoTensor& RankTwoTensor::operator=(const double & a)
{
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=a;
        }
    }
    return (*this);
}
RankTwoTensor& RankTwoTensor::operator=(const RankTwoTensor & a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: = applied to two diff rank2 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=a(i,j);
        }
    }
    return (*this);
}

//***************************
RankTwoTensor RankTwoTensor::operator+(const double & a) const
{
    RankTwoTensor temp(nDim,0.0);
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            temp(i,j)=(*this)(i,j)+a;
        }
    }
    return temp;
}
RankTwoTensor RankTwoTensor::operator+(const RankTwoTensor & a) const
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank2 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankTwoTensor temp(nDim,0.0);
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            temp(i,j)=(*this)(i,j)+a(i,j);
        }
    }

    return temp;
}
//***************
RankTwoTensor& RankTwoTensor::operator+=(const double & a)
{
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=(*this)(i,j)+a;
        }
    }
    return (*this);
}
RankTwoTensor& RankTwoTensor::operator+=(const RankTwoTensor & a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: += applied to two diff rank2 tensor\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=(*this)(i,j)+a(i,j);
        }
    }
    return (*this);
}
//************
RankTwoTensor RankTwoTensor::operator-(const double & a) const
{
    RankTwoTensor temp(nDim,0.0);
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            temp(i,j)=(*this)(i,j)-a;
        }
    }
    return temp;
}
RankTwoTensor RankTwoTensor::operator-(const RankTwoTensor & a) const
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: - applied to two diff rank2 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankTwoTensor temp(nDim,0.0);
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            temp(i,j)=(*this)(i,j)-a(i,j);
        }
    }

    return temp;
}
//***********************
RankTwoTensor& RankTwoTensor::operator-=(const double & a)
{
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=(*this)(i,j)-a;
        }
    }
    return (*this);
}
RankTwoTensor& RankTwoTensor::operator-=(const RankTwoTensor & a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: -= applied to two diff rank2 tensor\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=(*this)(i,j)-a(i,j);
        }
    }
    return (*this);
}

//*************************************
RankTwoTensor  RankTwoTensor::operator*(const double & a) const
{
    RankTwoTensor temp(nDim,0.0);
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            temp(i,j)=(*this)(i,j)*a;
        }
    }
    return temp;
}
vector<double> RankTwoTensor::operator*(const vector<double> &a) const
{
    if(nDim!=a.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: * applied to two diff rank12 tensor\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    vector<double> temp(nDim);
    for(int i=1;i<=nDim;i++)
    {
        temp[i-1]=0.0;
        for(int j=1;j<=nDim;j++)
        {
            temp[i-1]+=(*this)(i,j)*a[j-1];
        }
    }

    return temp;
}
RankTwoTensor RankTwoTensor::operator*(const RankTwoTensor & a) const
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: - applied to two diff rank2 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankTwoTensor temp(nDim,0.0);
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            temp(i,j)=0.0;
            for(int k=1;k<=nDim;k++)
            {
                temp(i,j)+=(*this)(i,k)*a(k,j);
            }
        }
    }
    return temp;
}
//****************************
RankTwoTensor & RankTwoTensor::operator*=(const double & a)
{
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=(*this)(i,j)*a;
        }
    }
    return (*this);
}
RankTwoTensor& RankTwoTensor::operator*=(const RankTwoTensor & a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: -= applied to two diff rank2 tensor\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    (*this)=(*this)*a;
    return (*this);
}

