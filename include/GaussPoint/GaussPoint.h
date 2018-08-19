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
// Created by walkandthinker on 19.08.18.
// Define gauss integration points and weights in AsFem

#ifndef ASFEM_GAUSSPOINT_H
#define ASFEM_GAUSSPOINT_H

#include "petsc.h"

void Int1D(int ngp,int &Lint,double (&gs)[125][4]);
void Int2D(int ngp,int &Lint,double (&gs)[125][4]);
void Int3D(int ngp,int &Lint,double (&gs)[125][4]);

#endif //ASFEM_GAUSSPOINT_H
