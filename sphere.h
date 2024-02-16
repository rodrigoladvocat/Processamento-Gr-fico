#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include "vec3.h"

class sphere
{
    public:
        sphere(const vec3& center, const double radius,
               const vec3& coefic_dif, const vec3& coefic_esp,
               const vec3& coefic_amb, const vec3& coefic_ref,
               const vec3& coefic_tr, const vec3& coefic_rug
               ): cent(center), rad(radius), coef_dif(coefic_dif),
               coef_esp(coefic_esp), coef_amb(coefic_amb), coef_ref(coefic_ref),
               coef_tr(coefic_tr), coef_rug(coefic_rug) {}

        vec3 cent;
        double rad;
        vec3 coef_dif;
        vec3 coef_esp;
        vec3 coef_amb;
        vec3 coef_ref;
        vec3 coef_tr;
        vec3 coef_rug;
        
};


#endif