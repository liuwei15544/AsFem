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
// define DofHanlder in AsFem

#include "DofHandler/DofHandler.h"

DofHandler::DofHandler()
{
    HasDofMap=false;
    HasBCDofMap=false;
    nDims=-1;nDofs=-1;nNodes=-1;nElmts=-1;
    nDofsPerNode=-1;nNodesPerElmts=-1;
    nDofsPerElmt=-1;nDofsPerBCElmt=-1;

    GlobalDofMap.clear();
    GlobalBCDofMap.clear();
    LeftBCDofs.clear();RightBCDofs.clear();
    BottomBCDofs.clear();TopBCDofs.clear();
    BackBCDofs.clear();FrontBCDofs.clear();

    nDofsPerBCElmt=-1;nNodesPerBCElmt=-1;
}


