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
// Created by walkandthinker on 31.08.18.
// define the free energy materials for AllenCahn model in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::FreeEnergyForAllenCahn(const double &conc, double &f, double &dfdc, double &d2fdc2)
{
    // For free energy and its derivative
    f=(conc*conc-1)*(conc*conc-1);
    dfdc=conc*conc*conc-conc;
    d2fdc2=3*conc*conc-1;
}