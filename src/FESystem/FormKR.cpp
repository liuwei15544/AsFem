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
// Created by walkandthinker on 19.08.18.
// Loop all the elements, form K and R for FE system

#include "FESystem/FESystem.h"

void FESystem::FormKR(const int &iState, const double dt, const double t, const double (&ctan)[2],
                      Mesh &mesh,DofHandler &dofHandler,BCSystem &bcSystem,
                      const Vec &U, const Vec &V,
                      Mat &AMATRIX, Vec &RHS,Mat &Proj)
{
    int iInd,jInd,e,i,j,k;
    PetscScalar value;

    // we can get the correct value on the ghosted node!
    VecScatterCreateToAll(U,&scatteru,&Useq);
    VecScatterBegin(scatteru,U,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatteru,U,Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(V,&scatterv,&Vseq);
    VecScatterBegin(scatterv,V,Vseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatterv,V,Vseq,INSERT_VALUES,SCATTER_FORWARD);

    // Initializing all the vector and matrix
    ZeroMatAndVec(iState,AMATRIX,Proj,RHS);

    // For parallel
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    int rankne=mesh.GetElmtsNum()/size;
    int eStart=rank*rankne;
    int eEnd=(rank+1)*rankne;
    if(rank==size-1) eEnd=mesh.GetElmtsNum();

    for(e=eStart;e<eEnd;++e)
    {
        nNodesPerElmt=mesh.GetNodesNumPerElmt();
        dofHandler.GetLocalDofMap(e+1,nDofsPerElmt,elDofsConn);

        for(i=1;i<=nNodesPerElmt;++i)
        {
            j=mesh.GetIthConnJthIndex(e+1,i);
            elConn[i-1]=j-1;
            elCoords[i-1][1]=mesh.GetIthNodeJthCoord(j,1);
            elCoords[i-1][2]=mesh.GetIthNodeJthCoord(j,2);
            elCoords[i-1][3]=mesh.GetIthNodeJthCoord(j,3);
            for(k=1;k<=nDofsPerNode;k++)
            {
                iInd=(i-1)*nDofsPerNode+k-1;
                jInd=elDofsConn[iInd]-1;
                VecGetValues(Useq,1,&jInd,&elU[iInd][0]);

                VecGetValues(Vseq,1,&jInd,&elU[iInd][1]);
            }
        }



        //ElmtLib(e+1,t,dt,uelnum,istate,nDims,nNodesPerElmt,nDofsPerElmt,elCoords,elU,ctan,localK,localRHS,localHist,localProj);


        //AssembleLocalToGlobal(AMATRIX,RHS)
    }

    //VecAssemblyBegin(RHS);
    //VecAssemblyEnd(RHS);

    //MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    //MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

    // delete scatter
    VecScatterDestroy(&scatteru);
    VecScatterDestroy(&scatterv);
    VecDestroy(&Useq);
    VecDestroy(&Vseq);
}

