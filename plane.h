#ifndef PLANE_H
#define PLANE_H

#include <iostream>
#include "vec3.h"

class plane
{
    public:
        plane(vec3 col, const vec3& plane_point, const vec3& plane_vector,
              double coefic_dif, double coefic_esp,
              double coefic_amb, double coefic_ref,
              double coefic_tr, double coefic_rug, double refrac_index):
              color(col), pp(plane_point), pv(plane_vector), coef_dif(coefic_dif),
              coef_esp(coefic_esp), coef_amb(coefic_amb), coef_ref(coefic_ref),
              coef_tr(coefic_tr), coef_rug(coefic_rug), refr_index(refrac_index) {}

        vec3 color;
        vec3 pp;
        vec3 pv;
        double coef_dif;
        double coef_esp;
        double coef_amb;
        double coef_ref;
        double coef_tr;
        double coef_rug;
        double refr_index;
};

#endif