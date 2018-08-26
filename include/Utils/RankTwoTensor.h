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

#ifndef ASFEM_RANKTWOTENSOR_H
#define ASFEM_RANKTWOTENSOR_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "petsc.h"

using namespace std;


class RankTwoTensor
{
public:
    enum InitMethod
    {
        InitZero,
        InitIdentity
    };
    RankTwoTensor(int dim,double val=0);
    RankTwoTensor(int dim,InitMethod method=InitZero);
    RankTwoTensor(double (&gradUx)[3],double (&gradUy)[3]);
    RankTwoTensor(double (&gradUx)[3],double (&gradUy)[3],double (&gradUz)[3]);
    // from voigt notation
    RankTwoTensor(double &s11,double &s22,double &s12); // 2D case
    RankTwoTensor(double &s11,double &s22,double &s33,double &s23,double &s13,double &s12);// 3D case

    void ZeroEntities();
    void IdentityEntities();

    inline int GetDim() const { return nDim;}

    // operator overload
    double operator()(int i,int j) const;
    double& operator()(int i,int j);

    RankTwoTensor & operator=(const double & a);
    RankTwoTensor & operator=(const RankTwoTensor & a);

    RankTwoTensor operator+(const double & a) const;
    RankTwoTensor operator+(const RankTwoTensor & a) const;
    RankTwoTensor & operator+=(const double & a);
    RankTwoTensor & operator+=(const RankTwoTensor & a);

    RankTwoTensor operator-(const double & a) const;
    RankTwoTensor operator-(const RankTwoTensor & a) const;
    RankTwoTensor & operator-=(const double & a);
    RankTwoTensor & operator-=(const RankTwoTensor & a);

    RankTwoTensor  operator*(const double & a) const;
    vector<double> operator*(const vector<double> &a) const;
    RankTwoTensor  operator*(const RankTwoTensor & a) const;

    RankTwoTensor & operator*=(const double & a);
    RankTwoTensor & operator*=(const RankTwoTensor & a);

    // math operation for rank-4 tensor


private:
    int nDim;
    double elements[3*3]={0.0};

};


#endif //ASFEM_RANKTWOTENSOR_H
