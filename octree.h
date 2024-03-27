#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <vector>
#include "triangles.h"
#include "vec3.h"
#include "ray.h"

class Octree {

private:

    double xL, xR;
    double yL, yR;
    double zL, zR;
    Octree* child[8];

    std::vector<int> triangle_surf;

public:

    Octree createNode(Octree &current_node, int child_index) const {    // criando o octante
        double xl, xr, yl, yr, zl, zr;

        if((child_index&(1)) == 0) xl = current_node.getXL(), xr = current_node.getMX();
        else xl = current_node.getMX(), xr = current_node.getXR();   // => colocando na segunda metade do eixo x se o filho tem indice ímpar

        if((child_index&(2)) == 0) yl = current_node.getYL(), yr = current_node.getMY();
        else yl = current_node.getMY(), yr = current_node.getYR();   // colocando na segunda metade do eixo y se o filho tem segundo bit menos significativo sendo 1

        if((child_index&(4)) == 0) zl = current_node.getZL(), zr = current_node.getMZ();
        else zl = current_node.getMZ(), zr = current_node.getZR();  // - - - segunda metade do eixo z se o filho tem terceiro bit menos significativo sendo 1

        Octree oct = Octree(xl, xr, yl, yr, zl, zr);

        return oct;
    }

    void insert_surface(int t_index){
        triangle_surf.push_back(t_index);
    }

    void insert(const vec3 triangle, vec3* points_list, int t_index, Octree* &current_node, int level){
        int child_index = 0;

        vec3 minT = getMin(triangle, points_list);
        vec3 maxT = getMax(triangle, points_list);

        double mx = current_node->getMX();
        double my = current_node->getMY();
        double mz = current_node->getMZ();

        if (((minT.e[0] <= mx) && (maxT.e[0] >= mx)) ||
            ((minT.e[1] <= my) && (maxT.e[1] >= my)) ||
            ((minT.e[2] <= mz) && (maxT.e[2] >= mz)) || level == 4) {
                
            current_node->insert_surface(t_index);
            return;
        }
 
        if(minT.e[0] > mx) child_index += 1;  // indo para a segunda metade do cubo no eixo x  => separando indices pares de impares
        if(minT.e[1] > my) child_index += 2;  // indo para a segunda metade do cubo no eixo y  => separando indices usando o segundo bit menos significativo
        if(minT.e[2] > mz) child_index += 4;  // indo para a segunda metade do cubo no eixo z  => separando indices usando o terceiro bit menos significativo

        Octree* next_node = current_node->getChild(child_index);

        if(next_node == nullptr) {
            Octree oct = createNode(*current_node, child_index);
            next_node = new Octree(0, 0, 0, 0, 0, 0);
            *next_node = oct;
            current_node->setChild(child_index, next_node);
        }

        insert(triangle, points_list, t_index, next_node, level+1);
        return;
    }

    bool hit_octant(const ray &ray, Octree& current_node) const {
        double t0 = 0.0;
        double tf = 5000;

        double tl, tr;

        if(ray.direction().e[0] != 0){
            tl = (current_node.getXL() - ray.origin().e[0]) / ray.direction().e[0];
            tr = (current_node.getXR() - ray.origin().e[0]) / ray.direction().e[0];
            if(tl > tr) std::swap(tl, tr); // => colocando como minimo o bound atingido antes

            t0 = std::max(t0, tl);
            tf = std::min(tf, tr);
        }

        if(ray.direction().e[1] != 0){
            tl = (current_node.getYL() - ray.origin().e[1]) / ray.direction().e[1];
            tr = (current_node.getYR() - ray.origin().e[1]) / ray.direction().e[1];
            if(tl > tr) std::swap(tl, tr);

            t0 = std::max(t0, tl);
            tf = std::min(tf, tr);
        }

        if(ray.direction().e[2] != 0){
            tl = (current_node.getZL() - ray.origin().e[2]) / ray.direction().e[2];
            tr = (current_node.getZR() - ray.origin().e[2]) / ray.direction().e[2];
            if(tl > tr) std::swap(tl, tr);

            t0 = std::max(t0, tl);
            tf = std::min(tf, tr);
        }

        return (t0 < tf); // condiçao para haver interseçao com octante
    }

    void find(const ray &ray, Octree* current_node, std::vector<int> &possible_hits) const {
        if(current_node == nullptr) return;
        if(!hit_octant(ray, *current_node)) return;

        const std::vector<int> tmp = current_node->gettriangle_surf();
        for(int idx: tmp) possible_hits.push_back(idx);

        for(int i = 0;i < 8; i++) {
            find(ray, current_node->getChild(i), possible_hits);
        }

    }
    
    Octree(double xl, double xr, double yl, double yr, double zl, double zr){
        xL = xl;
        xR = xr;
        yL = yl;
        yR = yr;
        zL = zl;
        zR = zr;

        for(int i=0;i<8;i++) child[i] = nullptr;
    }

    double getXL() const {return xL;}
    double getXR() const {return xR;}
    double getYL() const {return yL;}
    double getYR() const {return yR;}
    double getZL() const {return zL;}
    double getZR() const {return zR;}
    double getMX() const {return (xL + xR) / 2.0;}
    double getMY() const {return (yL + yR) / 2.0;}
    double getMZ() const {return (zL + zR) / 2.0;}

    Octree* getChild(int child_index) const {
    return child[child_index];
}

    void setChild(int child_index, Octree* node){
        child[child_index] = node;
    }

    std::vector<int> gettriangle_surf() const {
        return triangle_surf;
    }
};

#endif