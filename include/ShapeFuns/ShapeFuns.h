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
// Created by walkandthinker on 15.08.18.
// Define the standard shape function used in AsFem

#ifndef ASFEM_SHAPEFUNS_H
#define ASFEM_SHAPEFUNS_H

#include "petsc.h"

void Shp1D(const int &ndim,const int &nnodes,const double &xi,const double (&Coords)[27][4],double (&shp)[27][4],double &DetJac);
void Shp2D(const int &ndim,const int &nnodes,const double &xi,const double &eta,const double (&Coords)[27][4],double (&shp)[27][4],double &DetJac);
void Shp3D(const int &ndim,const int &nnodes,const double &xi,const double &eta,const double &zeta,const double (&Coords)[27][4],double (&shp)[27][4],double &DetJac);



#endif //ASFEM_SHAPEFUNS_H
