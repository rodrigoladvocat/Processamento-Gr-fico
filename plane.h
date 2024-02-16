#ifndef PLANE_H
#define PLANE_H

#include <iostream>
#include "vec3.h"

class plane
{
    public:
        plane(const vec3& plane_point, const vec3& plane_vector,
              const vec3& coefic_dif = vec3(0, 0, 0), const vec3& coefic_esp = vec3(0, 0, 0),
              const vec3& coefic_amb = vec3(0, 0, 0), const vec3& coefic_ref = vec3(0, 0, 0),
              const vec3& coefic_tr = vec3(0, 0, 0), const vec3& coefic_rug = vec3(0, 0, 0)):
              pp(plane_point), pv(plane_vector), coef_dif(coefic_dif),
              coef_esp(coefic_esp), coef_amb(coefic_amb), coef_ref(coefic_ref),
              coef_tr(coefic_tr), coef_rug(coefic_rug) {}

        vec3 pp;
        vec3 pv;
        vec3 coef_dif;
        vec3 coef_esp;
        vec3 coef_amb;
        vec3 coef_ref;
        vec3 coef_tr;
        vec3 coef_rug;
}

#endif