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

bool NonlinearSolver::Solve(Mesh &mesh,
                            DofHandler &dofHandler,
                            BCSystem &bcSystem,
                            EquationSystem &equationSystem,
                            ElementSystem &elementSystem,
                            FE &fe,
                            LinearSolver &linearSolver,
                            FESystemInfo &feSystemInfo)
{
    switch (SolverType)
    {
        case newtonraphson:
            return NewtonRaphson(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver,feSystemInfo);
            break;
        case modifynewtonraphson:
            break;
        case netownwithlinesearch:
            return NRWithLineSearch(mesh,dofHandler,bcSystem,equationSystem,elementSystem,fe,linearSolver,feSystemInfo);
            break;
        case arclength:
            break;
        default:
            break;
    }
}


