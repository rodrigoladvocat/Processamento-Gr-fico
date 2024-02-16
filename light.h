#ifndef LIGHT_H
#define LIGHT_H

#include "vec3.h"

class lights
{
    public:

    int n_lights;
    vec3* lights;

    lights(int num_lights, const vec3* llights){
        n_lights = num_lights;
        lights = llights;
    }
}

#endif