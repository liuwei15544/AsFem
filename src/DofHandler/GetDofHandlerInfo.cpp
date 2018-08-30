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
// get dofhandler information in AsFem

#include "DofHandler/DofHandler.h"


void DofHandler::GetLocalDofMap(const int e, int &ndofsperelmt, int (&rInd)[270]) const
{
    ndofsperelmt=nDofsPerElmt;

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is invalid!!! ***\n",e);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    int i;
    for(i=1;i<=nDofsPerElmt;i++)
    {
        rInd[i-1]=GlobalDofMap[(e-1)*nDofsPerElmt+i-1]-1;
    }
}

//*****************************
//*** For boundary mesh
void DofHandler::GetLocalBCDofMap(string sidename, const int e, int &ndofsperbcelmt, int (&Ind)[270]) const
{
    int i;
    ndofsperbcelmt=nDofsPerBCElmt;

    if(nDims==1)
    {
        if(sidename=="left")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=LeftBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="right")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=RightBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else
        {
            // TODO: add gmsh physical group support
        }
    }
    else if(nDims==2)
    {
        if(sidename=="left")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=LeftBCDofs[(e-1)*nDofsPerBCElmt+i-1];

            }
        }
        else if(sidename=="right")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=RightBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="bottom")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=BottomBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="top")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=TopBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else
        {
            // TODO: add gmsh physical group support
        }
    }
    else if(nDims==3)
    {
        if(sidename=="left")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=LeftBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="right")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=RightBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="bottom")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=BottomBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="top")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=TopBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="back")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=BackBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else if(sidename=="front")
        {
            for(i=1;i<=nDofsPerBCElmt;i++)
            {
                Ind[i-1]=FrontBCDofs[(e-1)*nDofsPerBCElmt+i-1];
            }
        }
        else
        {
            // TODO: add gmsh physical group support
        }
    }
}

//********************************
int DofHandler::GetBCSideDofsNum(string sidename) const
{
    if(!HasBCDofMap)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't get side's dofs number!!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        bc dof maps isn't generated!!!  ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(nDims==1)
    {
        if(sidename=="left")
        {
            return LeftBCDofs.size();
        }
        else if(sidename=="right")
        {
            return RightBCDofs.size();
        }
        else
        {
            // TODO: Use gmsh msh file information

            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: wrong side name of 1D case!!!   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
    else if(nDims==2)
    {
        if(sidename=="left")
        {
            return LeftBCDofs.size();
        }
        else if(sidename=="right")
        {
            return RightBCDofs.size();
        }
        else if(sidename=="bottom")
        {
            return BottomBCDofs.size();
        }
        else if(sidename=="top")
        {
            return TopBCDofs.size();
        }
        else
        {
            // TODO: Use gmsh msh file information

            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: wrong side name of 1D case!!!   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }
    else if(nDims==3)
    {
        if(sidename=="left")
        {
            return LeftBCDofs.size();
        }
        else if(sidename=="right")
        {
            return RightBCDofs.size();
        }
        else if(sidename=="bottom")
        {
            return BottomBCDofs.size();
        }
        else if(sidename=="top")
        {
            return TopBCDofs.size();
        }
        else if(sidename=="back")
        {
            return BackBCDofs.size();
        }
        else if(sidename=="front")
        {
            return FrontBCDofs.size();
        }
        else
        {
            // TODO: Use gmsh msh file information

            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: wrong side name of 1D case!!!   ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }

    return 0;
}


