#include <iostream>

#include "petsc.h"


//**********************
#include "Mesh/Mesh.h"
#include "InputSystem/InputSystem.h"
#include "EquationSystem/EquationSystem.h"

using namespace std;
int main(int args,char *argv[])
{
    PetscErrorCode ierr;
    cout<<"args="<<args<<": argv[0]="<<argv[0]<<endl;

    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    Mesh mesh;
    EquationSystem equationSystem;

    InputSystem inputSystem(args,argv);
    inputSystem.ReadMeshBlock(mesh);
    inputSystem.ReadDofsName(equationSystem);

    mesh.PrintMeshInfo();
    equationSystem.PrintSolutionNameMap();


    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}