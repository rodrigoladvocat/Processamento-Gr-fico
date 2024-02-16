#ifndef LIGHT_H
#define LIGHT_H

#include "vec3.h"

class light
{
    public:
        
        light(const vec3& localizacao, const vec3& cor): loc(localizacao), color(cor) {}

        vec3 loc;
        vec3 color;
    
};

#endif