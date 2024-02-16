#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "vec3.h"
#include <iostream>

class triangles
{
    public:
        triangles(vec3 col, const int num_t, vec3* point_list, vec3* triangle_list,
                  double coefic_dif, double coefic_esp,
                  double coefic_amb, double coefic_ref,
                  double coefic_tr, double coefic_rug
                  ): 
                  color(col), n_t(num_t), p_list(point_list), t_list(triangle_list),
                  coef_dif(coefic_dif), coef_esp(coefic_esp), coef_amb(coefic_amb), 
                  coef_ref(coefic_ref), coef_tr(coefic_tr), coef_rug(coefic_rug) {}
        
        vec3 color;
        int n_t;
        vec3* p_list;
        vec3* t_list;
        double coef_dif;
        double coef_esp;
        double coef_amb;
        double coef_ref;
        double coef_tr;
        double coef_rug;
};

#endif