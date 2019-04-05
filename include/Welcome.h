#ifndef ASFEM_WELCOME_H
#define ASFEM_WELCOME_H

#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;

void Welcome(double version=0.1)
{
    cout<<"*********************************************************"<<endl;
    cout<<"*** Welcome to use AsFem                              ***"<<endl;
    cout<<"*** A Simple Finite Element Method program            ***"<<endl;
    printf("*** Version: %5.1f                                    ***\n",version);
    cout<<"*** Author : walkandthinker                           ***"<<endl;
    cout<<"*** Contact: walkandthinker@gmail.com                 ***"<<endl;
    cout<<"*********************************************************"<<endl;
}

#endif 

