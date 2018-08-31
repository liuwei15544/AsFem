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
//*** Fill method
void RankTwoTensor::FillFromDispGradient(double (&gradUx)[3],double (&gradUy)[3])
{
    if(nDim!=2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: fill a non-2d grad tensor!!!    ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];
}
void RankTwoTensor::FillFromDispGradient(double (&gradUx)[3],double (&gradUy)[3],double (&gradUz)[3])
{
    if(nDim!=3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: fill a non-3d grad tensor!!!    ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];(*this)(1,3)=gradUx[2];
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];(*this)(2,3)=gradUy[2];
    (*this)(3,1)=gradUz[0];(*this)(3,2)=gradUz[1];(*this)(3,3)=gradUz[2];
}
//**************
void RankTwoTensor::ZeroEntities()
{
    for(int i=0;i<81;i++) elements[i]=0.0;
}
void RankTwoTensor::IdentityEntities()
{
    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            (*this)(i,j)=1.0*(i==j);
        }
    }
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
//****************************************
bool RankTwoTensor::operator==(const RankTwoTensor &a) const
{
    static const double tol=1.0e-13;

    if(GetDim()!=a.GetDim())
    {
        return false;
    }

    bool IsMatch=true;
    for(int i=1;i<=GetDim();i++)
    {
        for(int j=1;j<=GetDim();j++)
        {
            if(fabs((*this)(i,j)-a(i,j))>=tol)
            {
                IsMatch=false;
                return IsMatch;
            }
        }
    }
    return IsMatch;
}

//****************************************
//*** Other mathematical operations    ***
//****************************************
RankTwoTensor RankTwoTensor::transpose() const
{
    RankTwoTensor temp(GetDim(),0.0);
    for(int i=1;i<=GetDim();i++)
    {
        for(int j=1;j<=GetDim();j++)
        {
            temp(i,j)=(*this)(j,i);
        }
    }
    return temp;
}
//*************
void RankTwoTensor::vectorOuterProduct(vector<double> v1, vector<double> v2)
{
    if(v1.size()!=v2.size())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: v1 is not same as v2!!!         ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(v1.size()!=GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: vector isn't match rank2 tensor!***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    for(int i=1;i<=GetDim();i++)
    {
        for(int j=1;j<=GetDim();j++)
        {
            (*this)(i,j)=v1[i-1]*v2[j-1];
        }
    }
}

//*******************
double RankTwoTensor::doubleDot(const RankTwoTensor &a) const
{
    // sum over a_ij*b_ij(i.e. strain energy calculation)
    if(GetDim()!=a.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: (:) rank2 tensor size isn't match!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    double sum=0.0;
    for(int i=1;i<GetDim();i++)
    {
        for(int j=1;j<=GetDim();j++)
        {
            sum+=(*this)(i,j)*a(i,j);
        }
    }

    return sum;
}
//*** for trace
double RankTwoTensor::trace() const
{
    if(GetDim()==2)
    {
        return (*this)(1,1)+(*this)(2,2);
    }
    else
    {
        return (*this)(1,1)+(*this)(2,2)+(*this)(3,3);
    }
}
//*** for the determinant
double RankTwoTensor::det() const
{
    double determinant=0.0;
    if(GetDim()==2)
    {
        determinant=(*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
    }
    else
    {
        determinant=(*this)(1,1)*(*this)(2,2)*(*this)(3,3)
                   +(*this)(1,2)*(*this)(2,3)*(*this)(3,1)
                   +(*this)(1,3)*(*this)(2,1)*(*this)(3,2)
                   -(*this)(1,3)*(*this)(2,2)*(*this)(3,1)
                   -(*this)(1,2)*(*this)(2,1)*(*this)(3,3)
                   -(*this)(1,1)*(*this)(2,3)*(*this)(3,2);
    }

    return determinant;
}

//*** for the inverse
RankTwoTensor RankTwoTensor::inverse() const
{
    RankTwoTensor inv(GetDim(),0.0);
    double J=det();
    if(fabs(J)<1.0e-13)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: inv failed for singular tensor! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    if(GetDim()==2)
    {
        inv(1,1)= (*this)(2,2)/J;
        inv(1,2)=-(*this)(1,2)/J;
        inv(2,1)=-(*this)(2,1)/J;
        inv(2,2)= (*this)(1,1)/J;
    }
    else
    {
        // taken from wiki:
        //   https://en.wikipedia.org/wiki/Invertible_matrix
        double A= (*this)(2,2)*(*this)(3,3)-(*this)(2,3)*(*this)(3,2);
        double D=-(*this)(1,2)*(*this)(3,3)+(*this)(1,3)*(*this)(3,2);
        double G= (*this)(1,2)*(*this)(2,3)-(*this)(1,3)*(*this)(2,2);
        inv(1,1)=A/J;inv(1,2)=D/J;inv(1,3)=G/J;

        double B=-(*this)(2,1)*(*this)(3,3)+(*this)(2,3)*(*this)(3,1);
        double E= (*this)(1,1)*(*this)(3,3)-(*this)(1,3)*(*this)(3,1);
        double H=-(*this)(1,1)*(*this)(2,3)+(*this)(1,3)*(*this)(2,1);
        inv(2,1)=B/J;inv(2,2)=E/J;inv(2,3)=H/J;

        double C= (*this)(2,1)*(*this)(3,2)-(*this)(2,2)*(*this)(3,1);
        double F=-(*this)(1,1)*(*this)(3,2)+(*this)(1,2)*(*this)(3,1);
        double I= (*this)(1,1)*(*this)(2,2)-(*this)(1,2)*(*this)(2,1);
        inv(3,1)=C/J;inv(3,2)=F/J;inv(3,3)=I/J;
    }

    return inv;
}

//******************************************
//*** Print out tensor info
//******************************************
void RankTwoTensor::PrintTensor() const
{
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
    for(int i=1;i<=nDim;i++)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** ");
        for(int j=1;j<=nDim;j++)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%14.6e ",(*this)(i,j));
        }
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***\n");
    }
}

//****************************************
//*** For rank-4 tensor                ***
//****************************************
RankFourTensor RankTwoTensor::OuterProduct(const RankTwoTensor &b) const
{
    int i,j,k,l;
    if(GetDim()!=b.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: (x) rank2 tensor size isn't match!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankFourTensor temp(nDim,0.0);
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    temp(i,j,k,l)=(*this)(i,j)*b(k,l);
                }
            }
        }
    }

    return temp;
}

//****** For AIkBJl
RankFourTensor RankTwoTensor::AIkBJl(const RankTwoTensor &b) const
{
    int i,j,k,l;
    if(GetDim()!=b.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: (x) rank2 tensor size isn't match!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankFourTensor temp(nDim,0.0);
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    temp(i,j,k,l)=(*this)(i,k)*b(j,l);
                }
            }
        }
    }

    return temp;
}
//**** For AIlBJk
RankFourTensor RankTwoTensor::AIlBJk(const RankTwoTensor &b) const
{
    int i,j,k,l;
    if(GetDim()!=b.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: (x) rank2 tensor size isn't match!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankFourTensor temp(nDim,0.0);
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    temp(i,j,k,l)=(*this)(i,l)*b(j,k);
                }
            }
        }
    }

    return temp;
}
//*** For AJkBIl
RankFourTensor RankTwoTensor::AJkBIl(const RankTwoTensor &b) const
{
    int i,j,k,l;
    if(GetDim()!=b.GetDim())
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: (x) rank2 tensor size isn't match!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    RankFourTensor temp(nDim,0.0);
    for(i=1;i<=nDim;i++)
    {
        for(j=1;j<=nDim;j++)
        {
            for(k=1;k<=nDim;k++)
            {
                for(l=1;l<=nDim;l++)
                {
                    temp(i,j,k,l)=(*this)(j,k)*b(i,l);
                }
            }
        }
    }

    return temp;
}

//*** For Otimes
RankFourTensor RankTwoTensor::Otimes(const RankTwoTensor &b) const
{
    //AOtimeB=0.5*(AikBjl+AilBjk)
    return ((*this).AIkBJl(b)+(*this).AIlBJk(b))*0.5;
}