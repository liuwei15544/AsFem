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

    mesh.SetDim(3);

    mesh.SetXmin(0.0);
    mesh.SetXmax(10.0);
    mesh.SetYmin(0.0);
    mesh.SetYmax(2.0);
    mesh.SetZmin(1.0);
    mesh.SetZmax(2.0);

    mesh.SetNx(10);
    mesh.SetNy(20);
    mesh.SetNz(10);

    mesh.SetMeshType("hex27");
    mesh.CreateMesh();
    mesh.SaveMeshToVTU();



    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}