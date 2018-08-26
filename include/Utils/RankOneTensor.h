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

#ifndef ASFEM_RANKONETENSOR_H
#define ASFEM_RANKONETENSOR_H

#include <iostream>
#include <iomanip>

#include "petsc.h"



class RankOneTensor
{
public:
    RankOneTensor(int ndim=3,double val=0.0);

    int GetSize() const { return nDim;}
    void SetDim(int dim);
    void ZeroEntities();

    double operator()(int i) const;
    double& operator()(int i);

    RankOneTensor& operator=(double val);
    RankOneTensor&operator=(RankOneTensor &otherRankOneTensor);

    RankOneTensor operator+(double val);
    RankOneTensor operator+(RankOneTensor &otherRankOneTensor);

    RankOneTensor operator-(double val);
    RankOneTensor operator-(RankOneTensor &otherRankOneTensor);

    RankOneTensor operator*(double val);
    double operator*(RankOneTensor &otherRankOneTensor);


    void PrintTensor() const;

private:
    int nDim;
    double elements[3];
};


#endif //ASFEM_RANKONETENSOR_H
