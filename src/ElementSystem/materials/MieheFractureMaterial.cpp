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
// Created by walkandthinker on 29.09.18.
// for phase field fracture model(based on borden's work)

#include "ElementSystem/ElementSystem.h"

void ElementSystem::MieheFractureMaterial(const int &nDim,
                                           const double &conc,
                                           const RankTwoTensor &grad,
                                           const double &hist_old,
                                           double &hist,
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

    const double lamda=Parameters[0];
    const double mu=Parameters[1];
    const double l=Parameters[2];
    const double gc=Parameters[3];



    RankTwoTensor I(0.0);
    I.IdentityEntities();


    if(strainMode==small)
    {
        strain=0.5*(grad+grad.transpose());
    }
    else
    {
        RankTwoTensor F(0.0),Ft(0.0);
        I.IdentityEntities();

        F=grad+I;
        Ft=F.transpose();

        strain=0.5*(Ft*F-I);
    }






    RankTwoTensor eigvec(0.0);// store eigen vector
    vector<double> eigval;    // store eigen value

    RankFourTensor ProjPos=strain.PositiveProjectionTensor(eigval,eigvec);
    RankFourTensor I4Sym(0.0),ProjNeg(0.0);
    I4Sym.IdentityRank4Entities();

    ProjNeg=I4Sym-ProjPos;

    vector<RankTwoTensor> etens;
    RankTwoTensor temp(0.0);


    int i,j,k;
    for(k=1;k<=3;k++)
    {
        for(i=1;i<=3;i++)
        {
            for(j=1;j<=3;j++)
            {
                temp(i,j)=eigvec(i,k)*eigvec(j,k);
            }
        }
        etens.push_back(temp);
    }

    // Split the positive part and negative part
    vector<double> epos(3),eneg(3);

    for(i=1;i<=3;i++)
    {
        epos[i-1]=0.5*(abs(eigval[i-1])+eigval[i-1]);
        eneg[i-1]=0.5*(abs(eigval[i-1])-eigval[i-1]);
    }

    // for the negative and positive trace
    double etr=0.0;
    for(i=1;i<=3;i++)
    {
        etr+=eigval[i-1];
    }

    double etrpos=0.5*(abs(etr)+etr);
    double etrneg=0.5*(abs(etr)-etr);

    // Now , we can calculate the stress+ and stress-
    RankTwoTensor StressPos(0.0),StressNeg(0.0);
    for(i=1;i<=3;i++)
    {
        StressPos+=etens[i-1]*(lamda*etrpos+2.0*mu*epos[i-1]);
        StressNeg+=etens[i-1]*(lamda*etrneg+2.0*mu*eneg[i-1]);
    }

    // squres for energy calculation
    double pval,nval;
    pval=0.0;nval=0.0;
    for(i=1;i<=3;i++)
    {
        pval+=epos[i]*epos[i];
        nval+=eneg[i]*eneg[i];
    }

    // For free energy
    double PsiPos=0.5*lamda*etrpos*etrpos+mu*pval;
    double PsiNeg=0.5*lamda*etrneg*etrneg+mu*nval;

    if(PsiPos>hist)
    {
        hist=PsiPos;
    }
    else
    {
        hist=hist_old;
    }

    // Now, we can calculate the stress and jacobian

    stress=StressPos*(1-conc)*(1-conc)-StressNeg;

    dstressdc=StressPos*2.0*(1-conc)*(-1.0);

    // for the jacobian
    // dstress/dstrain=dstress_pos/dstrain_pos*dstrain_pos/dstrain
    //                +dstress_neg/dstrain_neg*dstrain_neg/dstrain
    // remark: dstrain_pos/dstrain=P_pos
    //         dstrain_neg/dstrain=P_neg

    double etrpos_sign,etrneg_sign;

    if(etr<0.0)
    {
        etrpos_sign=etr;
    }
    else
    {
        etrpos_sign=0.0;
    }

    if(-etr<0.0)
    {
        etrneg_sign=-etr;
    }
    else
    {
        etrneg_sign=0.0;
    }

    Jacobian=(1-conc)*(1-conc)*(lamda*I.OuterProduct(I)*etr+2.0*mu*ProjPos)
            +(lamda*I.OuterProduct(I)*etrneg_sign+2*mu*ProjNeg);


}

