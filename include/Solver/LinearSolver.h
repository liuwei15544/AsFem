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
// define linear solve for Ax=F

#ifndef ASFEM_LINEARSOLVER_H
#define ASFEM_LINEARSOLVER_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

using namespace std;

class LinearSolver
{
public:
    LinearSolver();
    bool InitSolver();
    bool Solve(Mat &A,Vec &x,Vec &F);

    void SetKSPMaxIterations(int maxiters);
    void SetKSPAbsoluteError(PetscReal atol);
    void SetKSPRelativeError(PetscReal rtol);
    void SetKSPDError(PetscReal dtol);

    // TODO: Let user choose the precondition type, not in command-line way!!!
    void SetPreconditionType(string pretype);

    int GetKSPIterations() const;
    int GetKSPMaxIterations() const;
    double GetKSPAbsoluteError() const;
    double GetKSPRelativeError() const;
    double GetKSPDivergenceError() const;

    string GetKSPPreconditionType() const;

    void Release();

    void PrintSolverInfo() const;

private:
    bool CurrentSolveFinished=false;
    bool IsSolverInit;
    KSP ksp;
    PetscErrorCode ierr;


    PC pc;

    string PreconditionType;
    PetscInt MaxIterations;
    PetscInt MaxGMRSRestat;
    PetscReal AbsoluteError,RelativeError,DError;
};


#endif //ASFEM_LINEARSOLVER_H
