#include <iostream>

#include "petsc.h"


//**********************
#include "FESystem/FESystem.h"
#include "FE/FE.h"
using namespace std;
int main(int args,char *argv[])
{
    PetscErrorCode ierr;

    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    FESystem feSystem;
    feSystem.Init(args,argv);


    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}