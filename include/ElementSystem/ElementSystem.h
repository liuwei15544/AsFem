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
// define the uel system in AsFem

#ifndef ASFEM_ELEMENTSYSTEM_H
#define ASFEM_ELEMENTSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>


#include "petsc.h"

//**********************************
//*** For AsFem's own header file
//**********************************
#include "ElementSystem/KernelBlockInfo.h"
#include "Utils/RankFourTensor.h"

using namespace std;

class ElementSystem
{
public:
    ElementSystem();


    KernelBlockInfo kernelBlockInfo;

public:
    void Init();
    void RunElmtLib(const int &iState,const int (&IX)[27],
                    const int &nDim,const int &nNodes,const int &nDofs,
                    const double &dt,const double &t,const double (&ctan)[2],
                    const double (&Coords)[27][4],const double (&U)[270][2],
                    double (&K)[270*270],double (&rhs)[270],double (&proj)[27][12+1]);

private:
    void SetUelIndex();
    void SetUmatIndex();


private:
    enum UelList
    {
        poisson,
        diffusion,
        cahnhilliard,
        mechanics,
        thermalmechanics,
        uel1,
        uel2,
        uel3,
        uel4,
        uel5,
        uel6,
        uel7,
        uel8,
        uel9,
        uel10
    };

    enum UmatList
    {
        emptyumat,
        linearelastic,
        neohookean,
        freeenergy,
        conductivity,
        umat1,
        umat2,
        umat3,
        umat4,
        umat5
    };

private:
    //**********************************************************
    //*** For UEL subroutines                                ***
    //**********************************************************
    void Poisson(const int &iState,const int (&IX)[27],
                 const int &nDim,const int &nNodes,const int &nDofs,
                 const double &dt,const double &t,const double (&ctan)[2],
                 const double (&Coords)[27][4],const double (&U)[270][2],
                 double (&K)[270*270],double (&rhs)[270],double (&proj)[27][12+1]);

    void SolidMechanics(const int &iState,const int (&IX)[27],
                 const int &nDim,const int &nNodes,const int &nDofs,
                 const double &dt,const double &t,const double (&ctan)[2],
                 const double (&Coords)[27][4],const double (&U)[270][2],
                 double (&K)[270*270],double (&rhs)[270],double (&proj)[27][12+1]);
    double ElasticityTensorComponent(int i,int k,const int &nDim,
                                     const RankFourTensor &elasticity_tensor,
                                     const int &iInd,const int &jInd,
                                     const double (&shp)[27][4]);

    void CahnHilliard(const int &iState,const int (&IX)[27],
                 const int &nDim,const int &nNodes,const int &nDofs,
                 const double &dt,const double &t,const double (&ctan)[2],
                 const double (&Coords)[27][4],const double (&U)[270][2],
                 double (&K)[270*270],double (&rhs)[270],double (&proj)[27][12+1]);

    void elmt01(const int &iState,const int (&IX)[27],
                const int &nDim,const int &nNodes,const int &nDofs,
                const double &dt,const double &t,const double (&ctan)[2],
                const double (&Coords)[27][4],const double (&U)[270][2],
                double (&K)[270*270],double (&rhs)[270],double (&proj)[27][12+1]);

private:
    void Projection(const int &nNodes,const double &xsj,
                    double (&shp)[27][4],const double (&value)[12],
                    double (&proj)[27][12+1]);

private:
    //**************************************************
    //*** For Umat system                            ***
    //**************************************************
    void MechanicsMaterials(const int &nDim,const RankTwoTensor &grad,
                           RankTwoTensor &strain,
                           RankTwoTensor &stress,
                           RankFourTensor &Jacobian);

    void LinearElasticMaterial(const int &nDim,const RankTwoTensor &grad,
                               RankTwoTensor &strain,
                               RankTwoTensor &stress,
                               RankFourTensor &Jacobian);

    void NeoHookeanMaterial(const int &nDim,const RankTwoTensor &grad,
                            RankTwoTensor &strain,
                            RankTwoTensor &stress,
                            RankFourTensor &Jacobian);


    void FreeEnergyMaterials(const double &conc,double &f,double &dfdc,double &d2fdc2);
    void FreeEnergyForCahnHilliard(const double &conc,double &f,double &dfdc,double &d2fdc2);

private:
    int ActiveUelIndex=-10000;
    int ActiveUmatIndex=-10000;
    bool IsInit=false;
    vector<double> Parameters;

};


#endif //ASFEM_ELEMENTSYSTEM_H
