#include <iostream>

#include "petsc.h"


//**********************
#include "Mesh/Mesh.h"
#include "InputSystem/InputSystem.h"
#include "EquationSystem/EquationSystem.h"
#include "ElementSystem/KernelBlockInfo.h"
#include "BCSystem/BCSystem.h"
#include "BCSystem/BCBlockInfo.h"

using namespace std;
int main(int args,char *argv[])
{
    PetscErrorCode ierr;
    cout<<"args="<<args<<": argv[0]="<<argv[0]<<endl;

    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    Mesh mesh;
    EquationSystem equationSystem;
    KernelBlockInfo kernelBlockInfo;
    BCSystem bcSystem;
    vector<BCBlockInfo> bcBlockList;

    InputSystem inputSystem(args,argv);
    inputSystem.ReadMeshBlock(mesh);
    inputSystem.ReadDofsName(equationSystem);
    inputSystem.ReadKernelBlock(kernelBlockInfo);
    inputSystem.ReadBoundaryBlock(bcBlockList);

    mesh.PrintMeshInfo();
    equationSystem.PrintSolutionNameMap();
    kernelBlockInfo.PrintKernelBlockInfo();

    bcSystem.InitFromBCBlockList(bcBlockList);
    bcSystem.PrintBCInfo();


    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}