#include <iostream>

#include "petsc.h"


//**********************
#include "FESystem/FESystem.h"

#include "Utils/RankOneTensor.h"
#include "Utils/RankTwoTensor.h"


using namespace std;
int main(int args,char *argv[])
{
    PetscErrorCode ierr;

    ierr=PetscInitialize(&args,&argv,NULL,NULL);CHKERRQ(ierr);

    //FESystem feSystem;
    //feSystem.Init(args,argv);
    //feSystem.Run();

    RankTwoTensor a(2,0.0);

    a(1,1)=1.0;
    a(2,1)=21;
    a(1,2)=12;
    a(2,2)=20.0;
    a.PrintTensor();

    RankTwoTensor inva=a.inverse();

    RankTwoTensor aT=a.transpose();

    aT.PrintTensor();

    cout<<"det aT="<<aT.det()<<", det a="<<a.det()<<", tr(a)="<<a.trace()<<endl;

    inva.PrintTensor();

    RankTwoTensor I=inva*a;

    I.PrintTensor();





    ierr=PetscFinalize();CHKERRQ(ierr);

    return ierr;
}