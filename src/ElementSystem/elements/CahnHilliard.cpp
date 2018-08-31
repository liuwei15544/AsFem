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
// Created by walkandthinker on 25.08.18.
// uel for CahnHilliard equation

#include "ElementSystem/ElementSystem.h"

#include "ShapeFuns/ShapeFuns.h"
#include "GaussPoint/GaussPoint.h"

void ElementSystem::CahnHilliard(const int &iState, const int (&IX)[27], const int &nDim, const int &nNodes,
                            const int &nDofs, const double &dt, const double &t, const double (&ctan)[2],
                            const double (&Coords)[27][4], const double (&U)[270][2], double (&K)[270*270],
                            double (&rhs)[270], double (&proj)[27][13])
{
    int i,j,k;
    int ngp,Lint,gpInd;
    int iInd,jInd;
    double shp[27][4],gs[125][4],xsj,JxW;
    double gradc[4],conc,cdot;
    double mu,gradmu[4];
    double xi,eta,zeta;
    double value[12];// currently, only 12 variable is allowed to be projected!
    double f,dfdc,d2fdc2;
    double M=1.0,kappa=2.0e-2;

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

    M=Parameters[0];
    kappa=Parameters[1];



    if(nDim==1)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: cahnhilliard only support 2d and 3d case!!!\n");
        PetscFinalize();
        abort();
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
        if(nDim==2)
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
        gradc[1]=0.0;gradc[2]=0.0;gradc[3]=0.0;conc=0.0;cdot=0.0;
        gradmu[1]=0.0;gradmu[2]=0.0;gradmu[3]=0.0;mu=0.0;
        for(i=0;i<nNodes;i++)
        {
            conc+=shp[i][0]*U[2*i  ][0];
            cdot+=shp[i][0]*U[2*i  ][1];

            mu  +=shp[i][0]*U[2*i+1][0];
            for(k=1;k<=nDim;k++)
            {
                gradc[k] +=shp[i][k]*U[2*i  ][0];
                gradmu[k]+=shp[i][k]*U[2*i+1][0];
            }
        }

        FreeEnergyMaterials(conc,f,dfdc,d2fdc2);

        // Calculate local K and RHS
        if(iState%3==0)
        {
            // calculate residual
            for(iInd=0;iInd<nNodes;iInd++)
            {
                // R_c
                rhs[2*iInd  ]+=cdot*shp[iInd][0]*JxW;
                // R_mu
                rhs[2*iInd+1]+=mu*shp[iInd][0]*JxW-dfdc*shp[iInd][0]*JxW;
                for(k=1;k<=nDim;k++)
                {
                    // For gradient term
                    // R_c
                    rhs[2*iInd  ]+=M*gradmu[k]*shp[iInd][k]*JxW;
                    // R_mu
                    rhs[2*iInd+1]+=-kappa*gradc[k]*shp[iInd][k]*JxW;
                }

                if(iState==6)
                {
                    // calculate jacobian
                    for(jInd=0;jInd<nNodes;jInd++)
                    {
                        // Kc,cdot
                        K[(2*iInd  )*nDofs+2*jInd  ]+=-shp[jInd][0]*shp[iInd][0]*JxW*ctan[1];

                        // Kmu,c
                        K[(2*iInd+1)*nDofs+2*jInd  ]+=d2fdc2*shp[jInd][0]*shp[iInd][0]*JxW*ctan[0];
                        // Kmu,mu
                        K[(2*iInd+1)*nDofs+2*jInd+1]+=-shp[jInd][0]*shp[iInd][0]*JxW*ctan[0];
                        for(k=1;k<=nDim;k++)
                        {
                            // Kc,mu
                            K[(2*iInd  )*nDofs+2*jInd+1]+=-M*shp[jInd][k]*shp[iInd][k]*JxW*ctan[0];

                            // Kmu,c
                            K[(2*iInd+1)*nDofs+2*jInd  ]+=kappa*shp[jInd][k]*shp[iInd][k]*JxW*ctan[0];
                        }
                    }
                }
            }
        }
        else if(iState==8)
        {
            // do projection: project gauss point's quantities to nodal point
            value[1-1]=conc;
            value[2-1]=mu;
            value[3-1]=f;
            value[4-1]=dfdc;
            value[5-1]=d2fdc2;
            value[6-1]=gradc[1];
            value[7-1]=gradc[2];
            value[8-1]=gradc[3];

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

