#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>
#include "vec3.h"

class sphere
{
    public:
        sphere(const vec3& col, const vec3& center, const double radius,
               double coefic_dif, double coefic_esp,
               double coefic_amb, double coefic_ref,
               double coefic_tr, double coefic_rug, double refrac_index
               ): color(col), cent(center), rad(radius), coef_dif(coefic_dif),
               coef_esp(coefic_esp), coef_amb(coefic_amb), coef_ref(coefic_ref),
               coef_tr(coefic_tr), coef_rug(coefic_rug), refr_index(refrac_index){}

        vec3 color;
        vec3 cent;
        double rad;
        double coef_dif;
        double coef_esp;
        double coef_amb;
        double coef_ref;
        double coef_tr;
        double coef_rug;
        double refr_index;
        
};


#endif