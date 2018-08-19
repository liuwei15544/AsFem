//
// Created by Y. Bai on 19.08.18.
//

#ifndef ASFEM_ELEMENTSYSTEM_H
#define ASFEM_ELEMENTSYSTEM_H


class ElementSystem
{
public:
    ElementSystem();

    enum ElementList
    {
        Poisson,
        Diffusion,
        CahnHilliard,
        Mechanics,
        ThermalMechanics,
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


};


#endif //ASFEM_ELEMENTSYSTEM_H
