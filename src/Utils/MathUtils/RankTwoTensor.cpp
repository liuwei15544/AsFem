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

#include "Eigen/Eigen/Dense"

//****************************************************
//*** For constructor
//****************************************************
RankTwoTensor::RankTwoTensor(double val)
{
    // for general cases(i.e. plane strain problem, this must be 3x3 )

    for(int i=0;i<9;i++) elements[i]=val;
}
//***
RankTwoTensor::RankTwoTensor(RankTwoTensor::InitMethod method)
{
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
    (*this)(1,1)=s11;(*this)(1,2)=s12;(*this)(1,3)=0.0;
    (*this)(2,1)=s12;(*this)(2,2)=s22;(*this)(2,3)=0.0;
    (*this)(3,1)=0.0;(*this)(3,2)=0.0;(*this)(3,3)=0.0;

}
RankTwoTensor::RankTwoTensor(double &s11,double &s22,double &s33,double &s23,double &s13,double &s12)
{
    // For 3D case
    (*this)(1,1)=s11;(*this)(1,2)=s12;(*this)(1,3)=s13;
    (*this)(2,1)=s12;(*this)(2,2)=s22;(*this)(2,3)=s23;
    (*this)(3,1)=s13;(*this)(3,2)=s23;(*this)(3,3)=s33;
}
//*** Construct from displacement gradient(or other gradient, here mainly for the F tensor)
RankTwoTensor::RankTwoTensor(double (&gradUx)[3],double (&gradUy)[3])
{
    // For 2D displacement gradient
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];(*this)(1,3)=0.0;
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];(*this)(2,3)=0.0;
    (*this)(3,1)=      0.0;(*this)(3,2)=      0.0;(*this)(3,3)=0.0;
}
RankTwoTensor::RankTwoTensor(double (&gradUx)[3],double (&gradUy)[3],double (&gradUz)[3])
{
    // For 3D displacement gradient
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];(*this)(1,3)=gradUx[2];
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];(*this)(2,3)=gradUy[2];
    (*this)(3,1)=gradUz[0];(*this)(3,2)=gradUz[1];(*this)(3,3)=gradUz[2];
}
//*** Fill method
void RankTwoTensor::FillFromDispGradient(double (&gradUx)[3],double (&gradUy)[3])
{
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];(*this)(1,3)=0.0;
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];(*this)(2,3)=0.0;
    (*this)(3,1)=      0.0;(*this)(3,2)=      0.0;(*this)(3,3)=0.0;
}
void RankTwoTensor::FillFromDispGradient(double (&gradUx)[3],double (&gradUy)[3],double (&gradUz)[3])
{
    (*this)(1,1)=gradUx[0];(*this)(1,2)=gradUx[1];(*this)(1,3)=gradUx[2];
    (*this)(2,1)=gradUy[0];(*this)(2,2)=gradUy[1];(*this)(2,3)=gradUy[2];
    (*this)(3,1)=gradUz[0];(*this)(3,2)=gradUz[1];(*this)(3,3)=gradUz[2];
}
//**************
void RankTwoTensor::ZeroEntities()
{
    for(int i=0;i<9;i++) elements[i]=0.0;
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
    RankTwoTensor temp(0.0);
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

    RankTwoTensor temp(0.0);
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
    RankTwoTensor temp(0.0);
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

    RankTwoTensor temp(0.0);
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
    RankTwoTensor temp(0.0);
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

    RankTwoTensor temp(0.0);
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
    RankTwoTensor temp(0.0);
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
    RankTwoTensor inv(0.0);
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

    RankFourTensor temp(0.0);
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

    RankFourTensor temp(0.0);
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

    RankFourTensor temp(0.0);
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

    RankFourTensor temp(0.0);
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
    RankFourTensor temp(0.0);
    temp=0.5*((*this).AIkBJl(b)+(*this).AIlBJk(b));
    return temp;
}

//*** for lhs * operator
RankTwoTensor operator*(const double &lhs, const RankTwoTensor &a)
{
    RankTwoTensor temp(0.0);
    for(int i=0;i<9;i++) temp.elements[i]=lhs*a.elements[i];
    return temp;
}


//************************************************
//*** For eigen value and eigen vector
//************************************************
void RankTwoTensor::GetEigenValueAndVector(vector<double> &eigenvalue,
                                                 RankTwoTensor &eigenvector)
{
    // taken from:
    // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3%C3%973_matrices
    if(nDim==2)
    {
        //
        // For two dimension case
        // use : http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
        // A=[a b;c d]
        //
        // Now I use Eigen's built-in eigen value and eigen vector solution
        // TODO:If there is well-documented algorithm, move to that simple and clean code
        //      Now I have to use Eigen, which is not a good choice for me.
        //      Include Eigen is too heavy for AsFem!!!
        //      But good thing is Eigen can handle non-symmetric case, even through
        //       the strain should be symmetric, any how, it works!

        Eigen::Matrix2d A;
        A<<(*this)(1,1),(*this)(1,2),
           (*this)(2,1),(*this)(2,2);

        Eigen::EigenSolver<Eigen::Matrix2d> es(A);

        if(es.eigenvalues().size()<1)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find eigen value for rank-2 !\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }

        Eigen::Matrix2d EigVec=es.eigenvectors().real();
        eigenvalue.clear();
        for(unsigned int i=0;i<es.eigenvalues().size();i++)
        {
            eigenvalue.push_back(es.eigenvalues().real()[i]);
        }

        eigenvector(1,1)=EigVec(0,0);eigenvector(2,1)=EigVec(1,0);
        eigenvector(1,2)=EigVec(0,1);eigenvector(2,2)=EigVec(1,1);

    }
    else if(nDim==3)
    {
        eigenvalue.clear();

        Eigen::Matrix3d A;
        A<<(*this)(1,1),(*this)(1,2),(*this)(1,3),
           (*this)(2,1),(*this)(2,2),(*this)(2,3),
           (*this)(3,1),(*this)(3,2),(*this)(3,3);

        Eigen::EigenSolver<Eigen::Matrix3d> es(A);

        if(es.eigenvalues().size()<1)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find eigen value for rank-2 !\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }

        Eigen::Matrix3d EigVec=es.eigenvectors().real();
        eigenvalue.clear();
        for(unsigned int i=0;i<es.eigenvalues().size();i++)
        {
            eigenvalue.push_back(es.eigenvalues().real()[i]);
        }

        eigenvector(1,1)=EigVec(0,0);
        eigenvector(2,1)=EigVec(1,0);
        eigenvector(3,1)=EigVec(2,0);

        eigenvector(1,2)=EigVec(0,1);
        eigenvector(2,2)=EigVec(1,1);
        eigenvector(3,2)=EigVec(2,1);

        eigenvector(1,3)=EigVec(0,2);
        eigenvector(2,3)=EigVec(1,2);
        eigenvector(3,3)=EigVec(2,2);
    }
}

//********************************************************
RankFourTensor RankTwoTensor::PositiveProjectionTensor(vector<double> &eigenvalue,
                                                       RankTwoTensor &eigenvec)
{
    // Algorithm is taken from:
    // C. Miehe and M. Lambrecht, Commun. Numer. Meth. Engng 2001; 17:337~353
    // https://onlinelibrary.wiley.com/doi/epdf/10.1002/cnm.404
    GetEigenValueAndVector(eigenvalue,eigenvec);

    // C=F^T F=lambda_i M_i where M_i=n_i x n_i
    //         lambda_i-->eigenvalue  n_i--> the eigenvector

    double epos[nDim],d[nDim]; // d is defined in Eq.(17)
    int a,b,i,j;
    for(i=1;i<=nDim;i++)
    {
        // here only use the positive part, which is defined as eps_pos=[abs(eps)+eps]/2.0
        epos[i-1]=0.5*(abs(eigenvalue[i])+eigenvalue[i]);

        d[i-1]=1.0;
        if(eigenvalue[i-1]<0.0) d[i-1]=0.0;
    }

    RankFourTensor ProjPos(0.0);

    // calculate Ma defined in Eq.(9)-2
    // Ma=n_a x n_a
    RankTwoTensor Ma(0.0);
    for(a=1;a<=nDim;a++)
    {
        // for Ma=n_a x n_a
        for(i=1;i<=nDim;i++)
        {
            for(j=1;j<=nDim;j++)
            {
                Ma(i,j)=eigenvec(i,a)*eigenvec(j,a);
                // here the eigenvec is a rank-2 tensor(nxn matrix):
                // i-th col is the i-th eigen vector
            }
        }

        // Eq.(19), first term on the right side
        ProjPos+=d[a-1]*Ma.OuterProduct(Ma);
    }

    // Now we calculate the Gab and Gba
    // We need a new rank-2 tensor Mb(same defination as Ma)
    RankTwoTensor Mb(0.0);
    RankFourTensor Gab(nDim,0.0);
    RankFourTensor Gba(nDim,0.0);
    double theta_ab;// defined in Eq.(21)-1
    const double tol=1.0e-14;
    for(a=1;a<=nDim;a++)
    {
        for(b=1;b<=nDim;b++)
        {
            //*********************************
            //*** For Ma and Mb
            for(i=1;i<=nDim;i++)
            {
                for(j=1;j<=nDim;j++)
                {
                    Ma(i,j)=eigenvec(i,a)*eigenvec(j,a);
                    Mb(i,j)=eigenvec(i,b)*eigenvec(j,b);
                }
            }

            Gab=Ma.AIkBJl(Mb)+Ma.AIlBJk(Mb);// Eq.(12)
            Gba=Mb.AIkBJl(Ma)+Mb.AIlBJk(Ma);// change the order of Eq.(12)

            // since only positive term is involved
            // e_a=0.5*(abs(lambda_a)+lambda_a)
            // P is defined as: 2dE/dC in Eq.(8)-2
            //  but 2dM/dC=(Gab+Gba)/(lambda_a-lambda_b)
            // E(C)=sum(e_a*M_a)
            // P=2dE(C)/dC=2(dE(C)/dM)*(dM/dC)

            if(abs(eigenvalue[a-1]-eigenvalue[b-1])<tol)
            {
                //if limit lambda_a to lambda_b in Eq.(24)
                theta_ab=0.5*(d[a-1]+d[b-1])/2.0;
            }
            else
            {
                // ea(lambda_a)=(1/2)*(lambda_a-1)
                // m=2 for green strain case in Eq.(16)
                theta_ab=0.5*(epos[a-1]-epos[b-1])/(eigenvalue[a-1]-eigenvalue[b-1]);// Eq.(21)-1
            }

            ProjPos+=theta_ab*(Gab+Gba);
        }
    }

    return ProjPos;
}