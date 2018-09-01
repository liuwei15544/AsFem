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
// Created by walkandthinker on 30.08.18.
// define the neohookean hyperelastic material behavior for mechanical problems in AsFem

#include "ElementSystem/ElementSystem.h"

void ElementSystem::NeoHookeanMaterial(const int &nDim,
                                       const RankTwoTensor &grad,
                                       RankTwoTensor &strain,
                                       RankTwoTensor &stress,
                                       RankFourTensor &Jacobian)
{
    static const double E=Parameters[0];
    static const double nu=Parameters[1];
    const double Lambda=E*nu/(1.+nu)/(1.-2.*nu);
    const double Mu=0.5*E/(1.+nu);
    if(Parameters.size()<2)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: for mechanics, you need E,nu!!!  ***\n");
        PetscFinalize();
        abort();
    }

    RankTwoTensor F(nDim,0.0),Ft(nDim,0.0),I(nDim,0.0);
    RankTwoTensor C(nDim,0.0),Cinv(nDim,0.0);
    RankTwoTensor pk2(nDim,0.0);
    double J;

    I.IdentityEntities();// set as diagonal unit tensor


    F=grad+I;// F=\nabla U+I
    Ft=F.transpose();
    C=Ft*F;
    Cinv=C.inverse();

    J=F.det();

    pk2=0.5*Lambda*(J*J-1.0)*Cinv+Mu*(I-Cinv);

    stress=pk2;

    // here the jacobian is inexact!!!
    Jacobian=Lambda*J*J*Cinv.OuterProduct(Cinv)
            +(2.0*Mu-Lambda*(J*J-1.0))*Cinv.Otimes(Cinv);


}