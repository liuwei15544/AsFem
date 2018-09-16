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
// Created by walkandthinker on 15.09.18.
// for thermal-mechanics coupled constitutive law

#include "ElementSystem/ElementSystem.h"

void ElementSystem::ThermalElasticMaterial(const int &nDim,
                                           const double &conc,
                                           const RankTwoTensor &grad,
                                           RankTwoTensor &strain,
                                           RankTwoTensor &stress,
                                           RankTwoTensor &dstressdc,
                                           RankFourTensor &Jacobian)
{
    if(Parameters.size()<4)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: for thermal-mechanics, you need:***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        E, nu, D and omega!!!           ***\n");
        PetscFinalize();
        abort();
    }

    const double E=Parameters[0];
    const double nu=Parameters[1];
    const double Omega=Parameters[3];


    RankTwoTensor I(nDim,0.0);
    I.IdentityEntities();


    if(strainMode==small)
    {

        strain=0.5*(grad+grad.transpose())-(Omega*conc/nDim)*I;
    }
    else
    {
        RankTwoTensor F(nDim,0.0),Ft(nDim,0.0);
        I.IdentityEntities();

        F=grad+I;
        Ft=F.transpose();

        strain=0.5*(Ft*F-I)-(Omega*conc/nDim)*I;
    }






    Jacobian.FillFromEandNu(E,nu);

    for(int i=1;i<=nDim;i++)
    {
        for(int j=1;j<=nDim;j++)
        {
            for(int k=1;k<=nDim;k++)
            {
                for(int l=1;l<=nDim;l++)
                    cout<<Jacobian(i,j,k,l)<<" ";
            }
            cout<<endl;
        }
    }
    cout<<endl<<endl;


    stress=Jacobian*strain;

    dstressdc=(-Omega/nDim)*Jacobian*I;


}

