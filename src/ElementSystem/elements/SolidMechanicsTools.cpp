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
// tools for solid mechanics(get the component for K matrix from rank-4 tensor)

#include "ElementSystem/ElementSystem.h"

double ElementSystem::ElasticityTensorComponent(int i, int k, const int &nDim,
                                                const RankFourTensor &elasticity_tensor,
                                                const int &iInd, const int &jInd, const double (&shp)[27][4])
{
    // i=1 for Ux,i=2 for Uy, i=3 for Uz
    // i.e. (i,k)=(1,1)===> For Kux,ux
    //      (i,k)=(1,2)===> For Kux,uy...
    //      iInd===>for test function
    //      jInd===>for trial function
    // When residual=stress_ij*dtest/dxj
    // we can have jacob=dresidual/duk
    //                  =d(stress_ij*dtest/dxj)/duk
    // stress itself=(C_ijmn*dum/dxn)
    // so we can have jacob=d(C_ijmn*dum/dxn*dtest/dxj)/duk
    // i->1 for sigma_xj, then k=2 for d()/duy=Kux,uy
    // K^{IJ}_{uiuk}=C_ijkl*N^{I}_{,j}*N^{J}_{,l}
    double val=0.0;
    int j,l;
    for(j=1;j<=nDim;j++)
    {
        for(l=1;l<=nDim;l++)
        {
            //val+=elasticity_tensor(i,j,k,l)*shp[jInd][l]*shp[iInd][j];
            //val+=elasticity_tensor.voigt(i,j,k,l)*shp[jInd][l]*shp[iInd][j];
            val+=elasticity_tensor(i,j,k,l)*shp[jInd][l]*shp[iInd][j];
        }
    }
    return val;
}

