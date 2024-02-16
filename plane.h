#ifndef PLANE_H
#define PLANE_H

#include <iostream>
#include "vec3.h"

class plane
{
    public:
        plane(const vec3& plane_point, const vec3& plane_vector,
              int coefic_dif, int coefic_esp,
              int coefic_amb, int coefic_ref,
              int coefic_tr, int coefic_rug):
              pp(plane_point), pv(plane_vector), coef_dif(coefic_dif),
              coef_esp(coefic_esp), coef_amb(coefic_amb), coef_ref(coefic_ref),
              coef_tr(coefic_tr), coef_rug(coefic_rug) {}

        vec3 pp;
        vec3 pv;
        int coef_dif;
        int coef_esp;
        int coef_amb;
        int coef_ref;
        int coef_tr;
        int coef_rug;
};

#endif