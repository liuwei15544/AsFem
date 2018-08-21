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
//  solve for Ax=F equations

#include "Solver/NonlinearSolver.h"

void NonlinearSolver::Solve(Mesh &mesh,
                            DofHandler &dofHandler,
                            BCSystem &bcSystem,
                            EquationSystem &equationSystem,
                            ElementSystem &elementSystem,
                            FE &fe,
                            LinearSolver &linearSolver)
{
    switch (SolverType)
    {
        case newtonraphson:
            if(!NewtonRaphson(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver))
            {

            }
            break;
        case modifynewtonraphson:
            break;
        case arclength:
            break;
        default:
            break;
    }
}


