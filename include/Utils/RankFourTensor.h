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

#ifndef ASFEM_RANKFOURTENSOR_H
#define ASFEM_RANKFOURTENSOR_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "petsc.h"

//***********************
#include "RankTwoTensor.h"

using namespace std;

class RankTwoTensor;

class RankFourTensor
{
public:
    enum InitMethod
    {
        InitZero,
        InitIdentity,
        InitRank4Identity
    };
public:
    RankFourTensor(int dim,double val);
    RankFourTensor(int dim,InitMethod method=InitZero);
    void FillFromEandNu(double E,double nu);// Young's modulus and poisson ratio
    void FillFromKandG(double K,double G);  // Bulk modulus and shear modulus

    void ZeroEntities();
    void IdentityEntities();
    void IdentityRank4Entities();

    inline int GetDim() const { return nDim;}
    //**********************************
    //*** operator overload          ***
    //**********************************
    double operator()(int i,int j,int k,int l) const;
    double& operator()(int i,int j,int k,int l);
    double voigt(int i,int j,int k,int l) const;

    RankFourTensor & operator=(const double & a);
    RankFourTensor & operator=(const RankFourTensor & a);

    RankFourTensor operator+(const double & a) const;
    RankFourTensor operator+(const RankFourTensor & a) const;
    RankFourTensor & operator+=(const double & a);
    RankFourTensor & operator+=(const RankFourTensor & a);

    RankFourTensor operator-(const double & a) const;
    RankFourTensor operator-(const RankFourTensor & a) const;
    RankFourTensor & operator-=(const double & a);
    RankFourTensor & operator-=(const RankFourTensor & a);

    RankFourTensor  operator*(const double & a) const;
    RankTwoTensor   operator*(const RankTwoTensor & a) const;
    RankFourTensor  operator*(const RankFourTensor & a) const;
    friend RankFourTensor operator*(const double &lhs,const RankFourTensor &a);

    RankFourTensor& operator*=(const double & a);
    RankFourTensor& operator*=(const RankFourTensor &a);


    bool operator==(const RankFourTensor &a) const;

    void PrintTensor() const;

private:
    int nDim;
    int Dim2,Dim3,Dim4;
    double elements[3*3*3*3]={0.0};
};


#endif //ASFEM_RANKFOURTENSOR_H
