#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include "vec3.h"

class sphere
{
    public:
        sphere(const vec3& center, const double radius,
               int coefic_dif, int coefic_esp,
               int coefic_amb, int coefic_ref,
               int coefic_tr, int coefic_rug
               ): cent(center), rad(radius), coef_dif(coefic_dif),
               coef_esp(coefic_esp), coef_amb(coefic_amb), coef_ref(coefic_ref),
               coef_tr(coefic_tr), coef_rug(coefic_rug) {}

        vec3 cent;
        double rad;
        int coef_dif;
        int coef_esp;
        int coef_amb;
        int coef_ref;
        int coef_tr;
        int coef_rug;
        
};


#endif