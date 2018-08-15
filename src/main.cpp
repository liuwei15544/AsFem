#include <iostream>

#include "petsc.h"


//**********************
#include "Mesh/Mesh.h"

using namespace std;
int main(int args,char *argv[])
{
    PetscErrorCode ierr;
    cout<<"args="<<args<<": argv[0]="<<argv[0]<<endl;

    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    Mesh mesh;

    mesh.SetDim(1);

    mesh.SetXmin(0.0);
    mesh.SetXmax(1.0);



    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}