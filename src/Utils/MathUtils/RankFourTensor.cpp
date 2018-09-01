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
// Created by walkandthinker on 29.08.18.
// Define the rank-4 tensor(i.e. jac[=dstress/dstrain]) in AsFem

#include "Utils/RankFourTensor.h"

//**************************************************
//*** For constructor                            ***
//**************************************************
RankFourTensor::RankFourTensor(int dim, double val)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid in rank4 tensor*\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    nDim=dim;
    Dim2=dim*dim;
    Dim3=dim*dim;
    Dim4=dim*dim*dim*dim;
    for(int i=0;i<81;i++) elements[i]=val;
}

RankFourTensor::RankFourTensor(int dim, RankFourTensor::InitMethod method)
{
    if(dim<1||dim>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%3d is invalid in rank4 tensor*\n",dim);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    nDim=dim;
    Dim2=dim*dim;
    Dim3=dim*dim;
    Dim4=dim*dim*dim*dim;

    if(method==InitZero)
    {
        for(int i=0;i<81;i++) elements[i]=0.0;
    }
    else if(method==InitIdentity)
    {
        for(int i=0;i<81;i++) elements[i]=0.0;
        for(int i=1;i<=nDim;i++)
        {
            (*this)(i,i,i,i)=1.0;
        }
    }
    else if(method==InitRank4Identity)
    {
        int iInd=0;
        for(int i=1;i<=nDim;i++)
        {
            for(int j=1;j<=nDim;j++)
            {
                for(int k=1;k<=nDim;k++)
                {
                    for(int l=1;l<=nDim;l++)
                    {
                        if(i==k && j==l)
                        {
                            elements[iInd]=1.0;
                        }
                        else
                        {
                            elements[iInd]=0.0;
                        }
                        iInd+=1;
                    }
                }
            }
        }
    }
}

//*****************************************************
//*** Operator overload                             ***
//*****************************************************
double RankFourTensor::operator()(int i, int j, int k, int l) const
{
    static int iInd,jInd;

    if(i<1||i>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range in rank4  ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(j<1||j>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%3d is out of range in rank4  ***\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(k<1||k>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: k=%3d is out of range in rank4  ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(l<1||l>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: l=%3d is out of range in rank4  ***\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    iInd=(i-1)*nDim+j;
    jInd=(iInd-1)*nDim+k;
    return elements[(jInd-1)*nDim+l-1];
}
double& RankFourTensor::operator()(int i, int j, int k, int l)
{
    static int iInd,jInd;

    if(i<1||i>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%3d is out of range in rank4  ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(j<1||j>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%3d is out of range in rank4  ***\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(k<1||k>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: k=%3d is out of range in rank4  ***\n",i);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(l<1||l>nDim)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: l=%3d is out of range in rank4  ***\n",j);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    iInd=(i-1)*nDim+j;
    jInd=(iInd-1)*nDim+k;
    return elements[(jInd-1)*nDim+l-1];
}
//*** For = operator
RankFourTensor& RankFourTensor::operator=(const double &a)
{
    for(int i=0;i<81;i++) elements[i]=a;
    return (*this);
}
RankFourTensor& RankFourTensor::operator=(const RankFourTensor &a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: = applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=0;i<Dim4;i++)
    {
        elements[i]=a.elements[i];
    }
    return (*this);
}
//*** For == operator
bool RankFourTensor::operator==(const RankFourTensor &a) const
{
    static const double tol=1.0e-13;
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:== applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=0;i<Dim4;i++)
    {
        if(fabs(elements[i]-a.elements[i])>=tol)
        {
            return false;
        }
    }
    return true;
}
//*** For + operator
RankFourTensor RankFourTensor::operator+(const double &a) const
{
    RankFourTensor temp(nDim,0.0);
    for(int i=0;i<Dim4;i++) temp.elements[i]=elements[i]+a;
    return temp;
}
RankFourTensor RankFourTensor::operator+(const RankFourTensor &a) const
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankFourTensor temp(nDim,0.0);
    for(int i=0;i<Dim4;i++) temp.elements[i]=elements[i]+a.elements[i];
    return temp;
}
RankFourTensor& RankFourTensor::operator+=(const double &a)
{
    for(int i=0;i<Dim4;i++) elements[i]=elements[i]+a;
    return (*this);
}
RankFourTensor& RankFourTensor::operator+=(const RankFourTensor &a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=0;i<Dim4;i++) elements[i]=elements[i]+a.elements[i];
    return (*this);
}
//*** For - operator
RankFourTensor RankFourTensor::operator-(const double &a) const
{
    RankFourTensor temp(nDim,0.0);
    for(int i=0;i<Dim4;i++) temp.elements[i]=elements[i]-a;
    return temp;
}
RankFourTensor RankFourTensor::operator-(const RankFourTensor &a) const
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankFourTensor temp(nDim,0.0);
    for(int i=0;i<Dim4;i++) temp.elements[i]=elements[i]-a.elements[i];
    return temp;
}
RankFourTensor& RankFourTensor::operator-=(const double &a)
{
    for(int i=0;i<Dim4;i++) elements[i]=elements[i]-a;
    return (*this);
}
RankFourTensor& RankFourTensor::operator-=(const RankFourTensor &a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=0;i<Dim4;i++) elements[i]=elements[i]-a.elements[i];
    return (*this);
}
//***************************
//*** For * operator
//***************************
RankFourTensor RankFourTensor::operator*(const double &a) const
{
    RankFourTensor temp(nDim,0.0);
    for(int i=0;i<Dim4;i++) temp.elements[i]=elements[i]*a;
    return temp;
}

RankTwoTensor RankFourTensor::operator*(const RankTwoTensor &a) const
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    // B_ij=C_ijkl*a_kl
    RankTwoTensor temp(nDim,0.0);
    int i,j,k,l;
    double val;
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            val=0.0;
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    val+=(*this)(i,j,k,l)*a(k,l);
                }
            }
            temp(i,j)=val;
        }
    }

    return temp;
}

RankFourTensor RankFourTensor::operator*(const RankFourTensor &a) const
{

    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    // E_ijkl=C_ijmn*D_mnkl
    RankFourTensor temp(nDim,0.0);
    int i,j,k,l,m,n;
    double val;
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    val=0.0;
                    for(m=1;m<=nDim;m++)
                    {
                        for(n=1;n<=nDim;n++)
                        {
                            val+=(*this)(i,j,m,n)*a(m,n,k,l);
                        }
                    }
                    temp(i,j,k,l)=val;
                }
            }
        }
    }
    return temp;
}
RankFourTensor& RankFourTensor::operator*=(const double &a)
{
    for(int i=0;i<Dim4;i++) elements[i]=elements[i]*a;
    return (*this);
}

RankFourTensor& RankFourTensor::operator*=(const RankFourTensor &a)
{
    if(nDim!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: + applied to two diff rank4 tensor*\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    (*this)=(*this)*a;
    return (*this);
}

//***************************************
//*** Fill method                     ***
//***************************************
void RankFourTensor::FillFromEandNu(double E, double nu)
{
    // fill out the rank-4 tensor from Young's modulus and nu
    const double Lambda=E*nu/((1.0+nu)*(1.0-2.0*nu));
    const double G=E/(2.0*(1.0+nu));
    int i,j,k,l;
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    (*this)(i,j,k,l)=
                            Lambda*(i==j)*(k==l)
                            +G*(i==k)*(j==l)
                            +G*(i==l)*(j==k);
                }
            }
        }
    }
}
//***
void RankFourTensor::FillFromKandG(double K, double G)
{
    // fill out the rank-4 tensor from Bulk modulus K and shear modulus G
    const double Lambda=K-2.0*G/3.0;
    int i,j,k,l;
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    (*this)(i,j,k,l)=
                            Lambda*(i==j)*(k==l)
                            +G*(i==k)*(j==l)
                            +G*(i==l)*(j==k);
                }
            }
        }
    }
}

//*************************************
//*** utils functions               ***
//*************************************
void RankFourTensor::ZeroEntities()
{
    for(int i=0;i<=Dim4;i++) elements[i]=0.0;
}

void RankFourTensor::IdentityEntities()
{
    for(int i=0;i<81;i++) elements[i]=0.0;
    for(int i=1;i<=nDim;i++)
    {
        (*this)(i,i,i,i)=1.0;
    }
}
void RankFourTensor::IdentityRank4Entities()
{
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            for(int k=1;k<=nDim;k++)
            {
                for(int l=1;l<=nDim;l++)
                {
                    if(i==k && j==l)
                    {
                        (*this)(i,j,k,l)=1.0;
                    }
                    else
                    {
                        (*this)(i,j,k,l)=0.0;
                    }
                }
            }
        }
    }
}

//******************************
RankFourTensor operator*(const double &lhs, const RankFourTensor &a)
{
    RankFourTensor temp(a.GetDim(),0.0);
    for(int i=0;i<81;i++) temp.elements[i]=lhs*a.elements[i];
    return temp;
}