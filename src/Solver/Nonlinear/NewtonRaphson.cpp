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
//  solve for Ax=F equations

#include "Solver/NonlinearSolver.h"

bool NonlinearSolver::NewtonRaphson(Mesh &mesh,
                                    DofHandler &dofHandler,
                                    BCSystem &bcSystem,
                                    EquationSystem &equationSystem,
                                    ElementSystem &elementSystem,
                                    FE &fe,
                                    LinearSolver &linearSolver,
                                    FESystemInfo &feSystemInfo)
{
    iters=0;IsConvergent=false;


    while(iters<=MaxIters && !IsConvergent)
    {
        bcSystem.ApplyDirichletBC(mesh,dofHandler,
                                  feSystemInfo.current_time,feSystemInfo.current_dt,
                                  equationSystem.U);


        fe.FormKR(feSystemInfo.iState,
                  feSystemInfo.current_dt,feSystemInfo.current_time,
                  feSystemInfo.ctan,
                  mesh,dofHandler,
                  elementSystem,equationSystem.U,
                  equationSystem.V,
                  equationSystem.AMATRIX,
                  equationSystem.RHS,
                  equationSystem.Proj);


        bcSystem.ApplyNeumannBC(mesh,dofHandler,
                                feSystemInfo.current_time,feSystemInfo.current_dt,
                                equationSystem.AMATRIX,equationSystem.RHS);
        bcSystem.ApplyConstraint(mesh,dofHandler,equationSystem.AMATRIX,equationSystem.RHS);


        if(iters==0 && feSystemInfo.currentstep==1)
        {
            MatSetOption(equationSystem.AMATRIX,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
        }
        linearSolver.Solve(equationSystem.AMATRIX,equationSystem.dU,equationSystem.RHS);


        VecNorm(equationSystem.RHS,NORM_2,&Rnorm);
        VecNorm(equationSystem.dU,NORM_2,&dUnorm);
        EnergyNorm=Rnorm*dUnorm;

        if(iters==0)
        {
            Rnorm0=Rnorm;
            dUnorm0=dUnorm;
            EnergyNorm0=EnergyNorm;
        }




        if(PrintIterationInfo && iters>0)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** iter=%6d                            ***\n",iters);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | R0|=%11.5e, | R|=%11.5e  ***\n",Rnorm0,Rnorm);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   |dU0|=%11.5e, |dU|=%11.5e  ***\n",dUnorm0,dUnorm);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | E0|=%11.5e, | E|=%11.5e  ***\n",EnergyNorm0,EnergyNorm);
        }

        equationSystem.UpdateUplusdU();
        //VecWAXPY(Vec w,PetscScalar a,Vec x,Vec y); w = a ∗ x + y
        VecWAXPY(equationSystem.V,-1.0,equationSystem.U0,equationSystem.U);//V=-U0+U
        //VecScale(Vec x, PetscScalar a); x = a ∗ x
        VecScale(equationSystem.V,feSystemInfo.ctan[1]);//V=V*(1.0/dt)


        iters+=1;

        if(ConvergenceCheck())
        {
            IsConvergent=true;
            break;
        }


    }

    if(!PrintIterationInfo)
    {
        if(feSystemInfo.jobtype=="static")
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** iter=%6d                            ***\n",iters);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | R0|=%11.5e, | R|=%11.5e  ***\n",Rnorm0,Rnorm);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   |dU0|=%11.5e, |dU|=%11.5e  ***\n",dUnorm0,dUnorm);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | E0|=%11.5e, | E|=%11.5e  ***\n",EnergyNorm0,EnergyNorm);
        }
        else
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** step=%6d, iter=%6d               ***\n",feSystemInfo.currentstep,iters);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | R0|=%11.5e, | R|=%11.5e  ***\n",Rnorm0,Rnorm);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   |dU0|=%11.5e, |dU|=%11.5e  ***\n",dUnorm0,dUnorm);
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***   | E0|=%11.5e, | E|=%11.5e  ***\n",EnergyNorm0,EnergyNorm);
        }
    }



    return IsConvergent;
}


