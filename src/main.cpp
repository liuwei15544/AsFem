#include <iostream>

#include "petsc.h"


//**********************
#include "Mesh/Mesh.h"
#include "InputSystem/InputSystem.h"

using namespace std;
int main(int args,char *argv[])
{
    PetscErrorCode ierr;
    cout<<"args="<<args<<": argv[0]="<<argv[0]<<endl;

    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    Mesh mesh;

    InputSystem inputSystem(args,argv);
    inputSystem.ReadMeshBlock(mesh);

    mesh.PrintMeshDetailInfo();


    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}