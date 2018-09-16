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
// for thermal-mechanics coupled problem

#include "ElementSystem/ElementSystem.h"

#include "ShapeFuns/ShapeFuns.h"
#include "GaussPoint/GaussPoint.h"

void ElementSystem::ThermalMechanics(const int &iState, const int (&IX)[27], const int &nDim, const int &nNodes,
                                     const int &nDofs, const double &dt, const double &t, const double (&ctan)[2],
                                     const double (&Coords)[27][4], const double (&U)[270][2], double (&K)[270*270],
                                     double (&rhs)[270], double (&proj)[27][13])
{
    int i,j,k;
    int ngp,Lint,gpInd;
    int iInd,jInd;
    double shp[27][4],gs[125][4],xsj,JxW;
    double gradUx[3],gradUy[3],gradUz[3];
    double gradc[3],gradSigmaH[3];
    double cdot,conc;
    double D,Omega,prefactor;
    double xi,eta,zeta;
    double value[12];// currently, only 12 variable is allowed to be projected!

    RankTwoTensor stress(0.0),dstressdc(0.0),strain(0.0),grad(0.0);
    RankFourTensor Jacobian(0.0);
    RankTwoTensor I(0.0);
    I.IdentityEntities();

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
        gradUx[0]=0.0;gradUx[1]=0.0;gradUx[2]=0.0;
        gradUy[0]=0.0;gradUy[1]=0.0;gradUy[2]=0.0;
        gradUz[0]=0.0;gradUz[1]=0.0;gradUz[2]=0.0;
        cdot=0.0;conc=0.0;
        gradc[0]=0.0;gradc[1]=0.0;gradc[2]=0.0;
        gradSigmaH[0]=0.0;gradSigmaH[1]=0.0;gradSigmaH[2]=0.0;
        for(i=0;i<nNodes;i++)
        {

            if(nDim==2)
            {
                conc+=shp[i][0]*U[3*i+2][0];
                cdot+=shp[i][0]*U[3*i+2][1];
                for(k=1;k<=nDim;k++)
                {
                    gradUx[k-1]+=shp[i][k]*U[3*i  ][0];
                    gradUy[k-1]+=shp[i][k]*U[3*i+1][0];
                    gradc[k-1] +=shp[i][k]*U[3*i+2][0];
                }
            }
            else if(nDim==3)
            {
                for(k=1;k<=nDim;k++)
                {
                    gradUx[k-1]+=shp[i][k]*U[4*i  ][0];
                    gradUy[k-1]+=shp[i][k]*U[4*i+1][0];
                    gradUz[k-1]+=shp[i][k]*U[4*i+2][0];
                    gradc[k-1] +=shp[i][k]*U[4*i+3][0];
                }
            }
        }


        if(nDim==2)
        {
            grad.FillFromDispGradient(gradUx,gradUy);
        }
        else if(nDim==3)
        {
            grad.FillFromDispGradient(gradUx,gradUy,gradUz);
        }


        ThermalElasticMaterial(nDim,conc,grad,strain,stress,dstressdc,Jacobian);

        // SigmaH=Sigma_ij*delta_ij/nDim
        D=Parameters[2];
        Omega=Parameters[3];

        prefactor=-(Omega/(3*3))*(Jacobian*I).trace();

        gradSigmaH[0]=prefactor*gradc[0];
        gradSigmaH[1]=prefactor*gradc[1];
        gradSigmaH[2]=prefactor*gradc[2];
        //prefactor=0.0;




        // Calculate local K and RHS
        if(iState%3==0)
        {
            // calculate residual
            for(iInd=0;iInd<nNodes;iInd++)
            {
                // R_c
                rhs[3*iInd+2]+=cdot*shp[iInd][0]*JxW;
                // R_u part
                // -Bt*Sigma
                if(nDim==2)
                {
                    for(k=1;k<=nDim;k++)
                    {
                        // R_ux
                        rhs[3*iInd  ]+=-stress(1,k)*shp[iInd][k]*JxW;
                        // R_uy
                        rhs[3*iInd+1]+=-stress(2,k)*shp[iInd][k]*JxW;
                        // R_c
                        rhs[3*iInd+2]+=D*gradc[k-1]*shp[iInd][k]*JxW
                                      -D*conc*Omega*gradSigmaH[k-1]*shp[iInd][k]*JxW;
                    }
                }
                else if(nDim==3)
                {
                    for(k=1;k<=nDim;k++)
                    {
                        // R_u
                        rhs[4*iInd  ]+=-stress(1,k)*shp[iInd][k]*JxW;
                        rhs[4*iInd+1]+=-stress(2,k)*shp[iInd][k]*JxW;
                        rhs[4*iInd+2]+=-stress(3,k)*shp[iInd][k]*JxW;
                        // R_c
                        rhs[4*iInd+3]+=D*gradc[k-1]*shp[iInd][k]*JxW
                                       -D*conc*Omega*gradSigmaH[k-1]*shp[iInd][k]*JxW;
                    }
                }

                if(iState==6)
                {
                    // calculate jacobian
                    for(jInd=0;jInd<nNodes;jInd++)
                    {
                        if(nDim==2)
                        {
                            // Kux,ux
                            K[(3*iInd  )*nDofs+3*jInd  ]+=ElasticityTensorComponent(1,1,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kux,uy
                            K[(3*iInd  )*nDofs+3*jInd+1]+=ElasticityTensorComponent(1,2,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];



                            // Kuy,ux
                            K[(3*iInd+1)*nDofs+3*jInd  ]+=ElasticityTensorComponent(2,1,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kuy,uy
                            K[(3*iInd+1)*nDofs+3*jInd+1]+=ElasticityTensorComponent(2,2,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];

                            // Kc,cdot
                            K[(3*iInd+2)*nDofs+3*jInd+2]+=-shp[jInd][0]*shp[iInd][0]*JxW*ctan[1];
                            for(k=1;k<=nDim;k++)
                            {
                                // K_ux,c
                                K[(3*iInd  )*nDofs+3*jInd+2]+=shp[jInd][0]*dstressdc(1,k)*shp[iInd][k]*JxW*ctan[0];
                                // K_uy,c
                                K[(3*iInd+1)*nDofs+3*jInd+2]+=shp[jInd][0]*dstressdc(2,k)*shp[iInd][k]*JxW*ctan[0];
                                // Kc,c
                                K[(3*iInd+2)*nDofs+3*jInd+2]+=-D*shp[jInd][k]*shp[iInd][k]*JxW*ctan[0]
                                                             +D*Omega*shp[jInd][0]*gradSigmaH[k-1]*shp[iInd][k]*JxW*ctan[0]
                                                             +D*conc*Omega*prefactor*shp[jInd][k]*shp[iInd][k]*JxW*ctan[0];
                            }

                        }
                        else if(nDim==3)
                        {
                            // Kux,ux
                            K[(4*iInd  )*nDofs+4*jInd  ]+=ElasticityTensorComponent(1,1,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kux,uy
                            K[(4*iInd  )*nDofs+4*jInd+1]+=ElasticityTensorComponent(1,2,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kux,uz
                            K[(4*iInd  )*nDofs+4*jInd+2]+=ElasticityTensorComponent(1,3,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];

                            // Kuy,ux
                            K[(4*iInd+1)*nDofs+4*jInd  ]+=ElasticityTensorComponent(2,1,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kuy,uy
                            K[(4*iInd+1)*nDofs+4*jInd+1]+=ElasticityTensorComponent(2,2,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kuy,uz
                            K[(4*iInd+1)*nDofs+4*jInd+2]+=ElasticityTensorComponent(2,3,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];

                            // Kuy,ux
                            K[(4*iInd+2)*nDofs+4*jInd  ]+=ElasticityTensorComponent(3,1,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kuy,uy
                            K[(4*iInd+2)*nDofs+4*jInd+1]+=ElasticityTensorComponent(3,2,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];
                            // Kuy,uz
                            K[(4*iInd+2)*nDofs+4*jInd+2]+=ElasticityTensorComponent(3,3,nDim,Jacobian,iInd,jInd,shp)*JxW*ctan[0];

                            // Kc,cdot
                            K[(4*iInd+3)*nDofs+4*jInd+3]+=-shp[jInd][0]*shp[iInd][0]*JxW*ctan[1];
                            for(k=1;k<=nDim;k++)
                            {
                                // K_ux,c
                                K[(4*iInd  )*nDofs+4*jInd+3]+=shp[jInd][0]*dstressdc(1,k)*shp[iInd][k]*JxW*ctan[0];
                                // K_uy,c
                                K[(4*iInd+1)*nDofs+4*jInd+3]+=shp[jInd][0]*dstressdc(2,k)*shp[iInd][k]*JxW*ctan[0];
                                // K_uz,c
                                K[(4*iInd+2)*nDofs+4*jInd+3]+=shp[jInd][0]*dstressdc(3,k)*shp[iInd][k]*JxW*ctan[0];
                                // Kc,c
                                K[(4*iInd+3)*nDofs+4*jInd+3]+=-D*shp[jInd][k]*shp[iInd][k]*JxW*ctan[0]
                                                              +D*Omega*shp[jInd][0]*gradSigmaH[k-1]*shp[iInd][k]*JxW*ctan[0]
                                                              +D*conc*Omega*prefactor*shp[jInd][k]*shp[iInd][k]*JxW*ctan[0];
                            }
                        }
                    }
                }
            }
        }
        else if(iState==8)
        {
            // do projection: project gauss point's quantities to nodal point
            if(nDim==2)
            {
                value[1-1]=stress(1,1);
                value[2-1]=stress(2,2);
                value[3-1]=stress(1,2);
                value[4-1]=strain(1,1);
                value[5-1]=strain(2,2);
                value[6-1]=strain(1,2);
            }
            else if(nDim==3)
            {
                value[1-1]=stress(1,1);
                value[2-1]=stress(2,2);
                value[3-1]=stress(3,3);
                value[4-1]=stress(2,3);
                value[5-1]=stress(1,3);
                value[6-1]=stress(1,2);

                value[ 7-1]=strain(1,1);
                value[ 8-1]=strain(2,2);
                value[ 9-1]=strain(3,3);
                value[10-1]=strain(2,3);
                value[11-1]=strain(1,3);
                value[12-1]=strain(1,2);
            }

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
