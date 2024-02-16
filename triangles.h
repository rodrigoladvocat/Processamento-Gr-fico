#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "vec3.h"
#include <iostream>

class triangles
{
    public:
        triangles(const int num_t, const vec3* point_list, const vec3* triangle_list,
                  const vec3& coefic_dif, const vec3& coefic_esp,
                  const vec3& coefic_amb, const vec3& coefic_ref,
                  const vec3& coefic_tr, const vec3& coefic_rug
                  ): 
                  n_t(num_t), p_list(point_list), t_list(triangle_list),
                  coef_dif(coefic_dif), coef_esp(coefic_esp), coef_amb(coefic_amb), 
                  coef_ref(coefic_ref), coef_tr(coefic_tr), coef_rug(coefic_rug) {}
        int n_t;
        vec3* p_list;
        vec3* t_list;
        vec3 coef_dif;
        vec3 coef_esp;
        vec3 coef_amb;
        vec3 coef_ref;
        vec3 coef_tr;
        vec3 coef_rug;
}

#endif