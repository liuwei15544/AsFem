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
// define nonlinear solve for Ax=F

#ifndef ASFEM_NONLINEARSOLVER_H
#define ASFEM_NONLINEARSOLVER_H

#include <iostream>
#include <string>

#include "petsc.h"

// For AsFem's own header file
#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "FE/FE.h"
#include "EquationSystem/EquationSystem.h"
#include "ElementSystem/ElementSystem.h"
#include "BCSystem/BCSystem.h"
#include "BCSystem/BCInfo.h"
#include "Solver/LinearSolver.h"
#include "FESystem/FESystemInfo.h"


using namespace std;

class NonlinearSolver
{
public:
    NonlinearSolver();


    void Init();
    void SetSolverType(string solvertype="newtonraphson");

    void SetPrintFlag(bool flag) {PrintIterationInfo=flag;}

    void SetMaxIters(int maxiters);
    int GetMaxIters() const { return MaxIters;}
    int GetCurrentIters() const { return iters;}

    void SetAtolOfResidual(double atol);
    void SetRtolOfResidual(double rtol);
    double GetAtolOfResidual() const { return Atol_R;}
    double GetRtolOfResidual() const { return Rtol_R;}

    void SetAtolOfdU(double atol);
    void SetRtolOfdU(double rtol);
    double GetAtolOfdU() const { return Atol_dU;}
    double GetRtolOfdU() const { return Rtol_dU;}

    double GetR0Norm() const { return Rnorm0;}
    double GetRnorm() const { return Rnorm;}

    double GetdU0Norm() const { return dUnorm0;}
    double GetdUNorm() const { return dUnorm;}

    double GetE0Norm() const { return EnergyNorm0;}
    double GetENorm() const { return EnergyNorm;}

    void SetAtolOfEnergy(double atol);
    void SetRtolOfEnergy(double rtol);
    double GetAtolOfEnergy() const { return Atol_E;}
    double GetRtolOfEnergy() const { return Rtol_E;}

    bool Solve(Mesh &mesh,
               DofHandler &dofHandler,
               BCSystem &bcSystem,
               EquationSystem &equationSystem,
               ElementSystem &elementSystem,
               FE &fe,
               LinearSolver &linearSolver,
               FESystemInfo &feSystemInfo);

    void PrintNonlinearSolverInfo() const;

private:
    bool NewtonRaphson(Mesh &mesh,
                       DofHandler &dofHandler,
                       BCSystem &bcSystem,
                       EquationSystem &equationSystem,
                       ElementSystem &elementSystem,
                       FE &fe,
                       LinearSolver &linearSolver,
                       FESystemInfo &feSystemInfo);

    bool NRWithLineSearch(Mesh &mesh,
                          DofHandler &dofHandler,
                          BCSystem &bcSystem,
                          EquationSystem &equationSystem,
                          ElementSystem &elementSystem,
                          FE &fe,
                          LinearSolver &linearSolver,
                          FESystemInfo &feSystemInfo);

    bool ModifyNewtonRaphson(Mesh &mesh,DofHandler &dofHandler,BCSystem &bcSystem,EquationSystem &equationSystem,LinearSolver &linearSolver);

    // TODO: arc-length method for multiple load step
    bool ArcLength();

private:
    bool ConvergenceCheck();

private:
    enum solvertype
    {
        newtonraphson,
        modifynewtonraphson,
        netownwithlinesearch,
        arclength
    };

private:
    int MaxIters,iters;
    solvertype SolverType;
    bool IsInit=false;
    double Atol_R,Rtol_R;   // absolute and relative error of residual
    double Atol_dU,Rtol_dU; // error of delta U
    long double Atol_E,Rtol_E;   // error of energy(=R*dU)

    bool IsConvergent;
    double      Rnorm,     Rnorm0;
    double     dUnorm,    dUnorm0;
    double EnergyNorm,EnergyNorm0;
    bool PrintIterationInfo=false;
};


#endif //ASFEM_NONLINEARSOLVER_H
