#ifndef TRIANGLES_H
#define TRIANGLES_H

#include "vec3.h"
#include <iostream>

using namespace std;

class triangles
{
    public:
        triangles(vec3 col, const int num_t, vec3* point_list, vec3* triangle_list,
                  double coefic_dif, double coefic_esp,
                  double coefic_amb, double coefic_ref,
                  double coefic_tr, double coefic_rug, double refrac_index
                  ): 
                  color(col), n_t(num_t), p_list(point_list), t_list(triangle_list),
                  coef_dif(coefic_dif), coef_esp(coefic_esp), coef_amb(coefic_amb), 
                  coef_ref(coefic_ref), coef_tr(coefic_tr), coef_rug(coefic_rug), refr_index(refrac_index) {}
        
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
        double refr_index;
};

point3 getMin(vec3 vertexes, vec3* points_list){  // pegando os minimos de um triangulo
    return point3( min(min(points_list[int(vertexes.e[0])].e[0], points_list[int(vertexes.e[1])].e[0]), points_list[int(vertexes.e[2])].e[0]),
                   min(min(points_list[int(vertexes.e[0])].e[1], points_list[int(vertexes.e[1])].e[1]), points_list[int(vertexes.e[2])].e[1]),
                   min(min(points_list[int(vertexes.e[0])].e[2], points_list[int(vertexes.e[1])].e[2]), points_list[int(vertexes.e[2])].e[2]));
}

point3 getMax(vec3 vertexes, vec3* points_list){  // pegando os maximos de um triangulo
    return point3( max(max(points_list[int(vertexes.e[0])].e[0], points_list[int(vertexes.e[1])].e[0]), points_list[int(vertexes.e[2])].e[0]),
                   max(max(points_list[int(vertexes.e[0])].e[1], points_list[int(vertexes.e[1])].e[1]), points_list[int(vertexes.e[2])].e[1]),
                   max(max(points_list[int(vertexes.e[0])].e[2], points_list[int(vertexes.e[1])].e[2]), points_list[int(vertexes.e[2])].e[2]));
}

#endif