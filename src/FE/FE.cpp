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
// Created by walkandthinker on 18.08.18.
// define the FE in AsFem

#include "FE/FE.h"

FE::FE()
{
    IsInit=false;

    nDims=0;nNodesPerElmt=-1;nDofsPerElmt=-1;nDofsPerNode=-1;
    current_elmt_id=-1;
}

