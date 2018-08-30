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
// Created by walkandthinker on 26.08.18.
// user defined element-01(uel) in AsFem

#include "ElementSystem/ElementSystem.h"

#include "ShapeFuns/ShapeFuns.h"
#include "GaussPoint/GaussPoint.h"

void ElementSystem::elmt01(const int &iState, const int (&IX)[27],const int &nDim,const int &nNodes,
                           const int &nDofs, const double &dt, const double &t, const double (&ctan)[2],
                           const double (&Coords)[27][4],const double (&U)[270][2], double (&K)[270*270],
                           double (&rhs)[270], double (&proj)[27][13])
{
    //**************************************************
    //*** This is an example for 2-dofs coupled case
    //*** Ru=grad(u)*grad(Ni)+grad(v)*grad(u)*Ni
    //*** Rv=grad(v)*grad(Ii)
    //**************************************************
    int i,j,k;
    int ngp,Lint,gpInd;
    int iInd,jInd;
    double shp[27][4],gs[125][4],xsj,JxW;
    double gradu[4],u;
    double gradv[4],v;
    double xi,eta,zeta;
    double value[12];// currently, only 12 variable is allowed to be projected!

    //******************************
    //*** Initializing
    //******************************
    if(iState%3==0)
    {
        // initializing local rhs
        for(i=0;i<nDofs;i++)
        {
            rhs[i]=0.0;
            if(iState==6)
            {
                // initializing locak k
                for(j=0;j<nDofs;j++)
                {
                    K[i*nDofs+j]=0.0;
                }
            }
        }
    }
    else if(iState==8)
    {
        // initializing projection array
        for(i=0;i<27;i++)
        {
            if(i<12) value[i]=0.0;
            for(j=0;j<=12;j++)
                proj[i][j]=0.0;
        }
    }
    //******************************************************

    ngp=3;

    if(nDim==1)
    {
        Int1D(ngp,Lint,gs);
    }
    else if(nDim==2)
    {
        Int2D(ngp,Lint,gs);
    }
    else if(nDim==3)
    {
        Int3D(ngp,Lint,gs);
    }

    for(gpInd=0;gpInd<Lint;gpInd++)
    {
        if(nDim==1)
        {
            xi=gs[gpInd][1];
            Shp1D(nDim,nNodes,xi,Coords,shp,xsj);
        }
        else if(nDim==2)
        {
            xi=gs[gpInd][1];
            eta=gs[gpInd][2];
            Shp2D(nDim,nNodes,xi,eta,Coords,shp,xsj);
        }
        else if(nDim==3)
        {
            xi=gs[gpInd][1];
            eta=gs[gpInd][2];
            zeta=gs[gpInd][3];
            Shp3D(nDim,nNodes,xi,eta,zeta,Coords,shp,xsj);
        }

        JxW=gs[gpInd][0]*xsj;


        // Calculate physical quantities on each gauss point
        gradu[1]=0.0;gradu[2]=0.0;gradu[3]=0.0;u=0.0;
        gradv[1]=0.0;gradv[2]=0.0;gradv[3]=0.0;v=0.0;
        for(i=0;i<nNodes;i++)
        {
            u+=shp[i][0]*U[2*i  ][0];
            v+=shp[i][0]*U[2*i+1][0];
            for(k=1;k<=nDim;k++)
            {
                gradu[k]+=shp[i][k]*U[2*i  ][0];
                gradv[k]+=shp[i][k]*U[2*i+1][0];
            }
        }

        // Calculate local K and RHS
        if(iState%3==0)
        {
            // calculate residual
            for(iInd=0;iInd<nNodes;iInd++)
            {
                for(k=1;k<=nDim;k++)
                {
                    // Ru
                    rhs[2*iInd  ]+=gradu[k]*shp[iInd][k]*JxW
                                +(gradv[k]*gradu[k])*shp[iInd][0]*JxW;
                    // Rv
                    rhs[2*iInd+1]+=gradv[k]*shp[iInd][k]*JxW;
                }

                if(iState==6)
                {
                    // calculate jacobian
                    for(jInd=0;jInd<nNodes;jInd++)
                    {
                        for(k=1;k<=nDim;k++)
                        {
                            // Kuu
                            K[(2*iInd  )*nDofs+2*jInd  ]+=-shp[jInd][k]*shp[iInd][k]*JxW*ctan[0]
                                                   -(gradv[k]*shp[jInd][k])*shp[iInd][0]*ctan[0]*JxW;
                            // Kuv
                            K[(2*iInd  )*nDofs+2*jInd+1]+=-(shp[jInd][k]*gradu[k])*shp[iInd][0]*ctan[0]*JxW;

                            // Kvv
                            K[(2*iInd+1)*nDofs+2*jInd+1]+=-shp[jInd][k]*shp[iInd][k]*ctan[0]*JxW;
                        }
                    }
                }
            }
        }
        else if(iState==8)
        {
            // do projection: project gauss point's quantities to nodal point
            value[1-1]=u;
            value[2-1]=gradu[1];
            value[3-1]=gradu[2];
            value[4-1]=gradu[3];
            value[5-1]=v;
            value[6-1]=gradv[1];
            value[7-1]=gradv[2];
            value[8-1]=gradv[3];

            Projection(nNodes,xsj,shp,value,proj);
        }
        else if(iState==5)
        {
            // TODO: initializing history value
        }
        else if(iState==10)
        {
            // TODO: update history value
        }
    }
}
