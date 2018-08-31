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
// convergence criteria for nonlinear iteration

#include "Solver/NonlinearSolver.h"

bool NonlinearSolver::ConvergenceCheck()
{
    if(Rnorm<Atol_R || EnergyNorm<Atol_E)
    {
        return true;
    }
    return false;
}

