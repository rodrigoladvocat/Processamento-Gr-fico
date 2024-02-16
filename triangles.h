#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "vec3.h"
#include <iostream>

class triangles
{
    public:
        triangles(const int num_t, vec3* point_list, vec3* triangle_list,
                  int coefic_dif, int coefic_esp,
                  int coefic_amb, int coefic_ref,
                  int coefic_tr, int coefic_rug
                  ): 
                  n_t(num_t), p_list(point_list), t_list(triangle_list),
                  coef_dif(coefic_dif), coef_esp(coefic_esp), coef_amb(coefic_amb), 
                  coef_ref(coefic_ref), coef_tr(coefic_tr), coef_rug(coefic_rug) {}
        int n_t;
        vec3* p_list;
        vec3* t_list;
        int coef_dif;
        int coef_esp;
        int coef_amb;
        int coef_ref;
        int coef_tr;
        int coef_rug;
};

#endif