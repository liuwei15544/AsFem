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
// for the petsc based math utils in equation system of AsFem

#include "EquationSystem/EquationSystem.h"

//************************************************
double EquationSystem::GetVecNorm(Vec &v)
{
    double norm;
    VecNorm(v,NORM_2,&norm);
    return norm;
}
//****************************
double EquationSystem::GetMatNorm(Mat &a)
{
    double norm;
    MatNorm(a,NORM_FROBENIUS,&norm);
    return norm;
}
//***************************************************

//******************************
double EquationSystem::GetRNorm()
{
    return GetVecNorm(RHS);
}
//*****************************
double EquationSystem::GetdUNorm()
{
    return GetVecNorm(dU);
}
//*****************************
double EquationSystem::GetUNorm()
{
    return GetVecNorm(U);
}
//***************************
double EquationSystem::GetEnergyNorm()
{
    double normR=GetRNorm();
    double normdU=GetdUNorm();
    return normR*normdU;
}
//**************************
double EquationSystem::GetKNrom()
{
    return GetMatNorm(AMATRIX);
}