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
// define the free energy materials for CahnHilliard model in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::FreeEnergyForCahnHilliard(const double &conc, double &f, double &dfdc, double &d2fdc2)
{
    // For free energy and its derivative
    f=100*conc*conc*(1-conc)*(1-conc);
    dfdc=200*conc*(conc -1)*(2*conc-1); // chemical potential
    d2fdc2=1200*conc*conc-1200*conc+200;
}