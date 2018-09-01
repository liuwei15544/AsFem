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
    static const double Lambda=Parameters[0];
    static const double Mu=Parameters[1];
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

    pk2=I*Mu-Cinv*Mu+Cinv*Lambda*log(J);

    //stress=F*pk2;// use pk1 not pk2

    stress=pk2;

    // here the jacobian is inexact!!!
    Jacobian=(Cinv.Otimes(Cinv))*2.0*Mu
            +Cinv.OuterProduct(Cinv)*Lambda
            -Cinv.Otimes(Cinv)*2.0*Lambda*log(J);
}