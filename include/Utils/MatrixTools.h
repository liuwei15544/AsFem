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
// Define the math functions for matrix operation

#ifndef ASFEM_MATRIXTOOLS_H
#define ASFEM_MATRIXTOOLS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "petsc.h"

using namespace std;

double det(const double (&A)[2][2]);
double det(const double (&A)[3][3]);

void inv(const double (&A)[2][2],double (&XA)[2][2]);
void inv(const double (&A)[3][3],double (&XA)[3][3]);

#endif //ASFEM_MATRIXTOOLS_H
