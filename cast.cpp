#include <iostream>
#include "vec3.h"
#include "ray.h"
#include "color.h"
#include "sphere.h"
#include "plane.h"
#include "triangles.h"
#include "light.h"
#include <vector>
#include <cmath>

#define EPSILON 2.22045e-016  // menor diferença entre dois doubles

vec3 triangle_normal(const point3& vertices, point3 *points_list)
{
    // finding the plane that contains the triangle

    int first = vertices.x();
    int second = vertices.y();
    int third = vertices.z();

    point3 A = points_list[first];
    point3 B = points_list[second];
    point3 C = points_list[third];

    vec3 vec_1 = A - B;
    vec3 vec_2 = C - B;

    vec3 t_normal = cross(vec_1, vec_2);  // vetor normal ao plano que é descrito pelos três vertices do triangulo

    return unit_vector(t_normal);
}

vec3 sphere_normal(const point3& center, const point3& intersection){
    return unit_vector(intersection - center);
}

double hit_sphere(sphere sphere, const ray& r) {
    vec3 oc = r.origin() - sphere.cent;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - sphere.rad*sphere.rad;
    auto delta = b*b - 4*a*c;

    if (delta < 0) return -1;

    return (std::min((- b - delta)/(2*a), (- b + delta)/(2*a)));  // retornando o menor valor de t para r(t) = v*t + o
}

double hit_plane(const point3& plane_point, const vec3& plane_vector, const ray& r){

    double denom = dot(r.direction(), plane_vector);
    if (denom != 0)
    {
        double t = (dot(plane_vector, plane_point) - dot(r.origin(), plane_vector)) / denom;
        return t;  
    }
    return -1;

}

double hit_triangle(const point3& vertices, point3 *points_list, const ray& r){
    
    int first = vertices.x();
    int second = vertices.y();
    int third = vertices.z();

    point3 A = points_list[first];
    point3 B = points_list[second];
    point3 C = points_list[third];

    vec3 t_normal = triangle_normal(vertices, points_list);  // vetor normal ao plano que é descrito pelos três vertices do triangulo

    double t = hit_plane(A, t_normal, r);

    if (t < 1) return -1; // o raio nao tem interseção com o plano

    point3 P = r.at(t);  // ponto de interseçao com o plano = > P = t*r + origem

    // para as areas dos triangulos formados pelas combinaçoes de vertices
    double aABC = 0.5 * cross((A - B), (C - B)).length();
    double aPBC = 0.5 * cross((B - P), (C - P)).length();
    double aPBA = 0.5 * cross((B - P), (A - P)).length();
    double aPAC = 0.5 * cross((A - P), (C - P)).length();

    // aPBC = aPBC / aABC;
    // aPBA = aPBA / aABC;
    // aPAC = aPAC / aABC;

    double sum = aPBC + aPBA + aPAC;

    if (fabs(aABC - sum) <= EPSILON) return t;
    else return -1;

}

// formula para o vetor normal à superficie e vetor incidente normalizados
// vetor incidente apontando da interseçao até a fonte de luz
vec3 reflected_light(vec3 normal_vector, vec3 incident_vector){
    vec3 product = 2 * dot(unit_vector(normal_vector), unit_vector(incident_vector)) * unit_vector(normal_vector);

    return unit_vector(product - incident_vector);
}

vec3 phong_equation(double ka, const color& ia, int n_lights, std::vector<light> light, 
                    point3 intersection, color od, double kd, vec3 normal_vector,
                    double ks, point3 camera, double krug){
    vec3 ambient = ka * ia;

    // somatorio:

    vec3 sum = ambient;
    for (int i = 0; i < n_lights; i++){

        vec3 li = unit_vector(light[i].loc - intersection); // intersecao -> fonte de luz
        vec3 reflection = reflected_light(normal_vector, li); // calculando vetor da luz refletida na sup.
        vec3 v = unit_vector(camera - intersection); // intersecao -> observador (camera)

        color iLn = light[i].color;
        
        vec3 first_factor = iLn * od * kd * dot(normal_vector, li);

        vec3 second_factor = iLn * ks * pow(dot(reflection, v), krug);
        
        sum = sum + first_factor + second_factor;
    }

    return sum;
}

color ray_color(const ray& r, triangles malha, std::vector<light> l_list, color filter, point3 camera) {
    double min_t = -1;
    vec3 min_normal;
    color min_t_color = color(0, 0, 0);

    double cd = 0, ce = 0, ca = 0, cr = 0, ct = 0, crug = 0;

   // declarando cada objeto
   
    sphere s1 = sphere(color(0, 0, 0), point3(2,0,0), 0.2, 
    0, 0, 0, 0, 0, 0);  
    sphere s2 = sphere(color(0, 0, 0), point3(5,5,0), 0.4, 
    0, 0, 0, 0, 0, 0);

    plane p1 = plane(color(0, 0, 0), point3(7, 2, 2), vec3(1, 0.2, 0.3),
    0, 0, 0, 0, 0, 0);

    std::vector<sphere> s_list;
    std::vector<plane> p_list;
    
    s_list.push_back(s1);
    s_list.push_back(s2);

    p_list.push_back(p1);

    int n_spheres = s_list.size();
    int n_planes = p_list.size();

    double t;
    for(int i = 0; i < malha.n_t; i++){
        t = hit_triangle(malha.t_list[i], malha.p_list, r);
        if (t >= 1){
            if (t < min_t || min_t == -1){
                min_t = t;
                min_normal = triangle_normal(malha.t_list[i], malha.p_list);
                cd = malha.coef_dif;
                ce = malha.coef_esp;
                ca = malha.coef_amb;
                cr = malha.coef_ref;
                ct = malha.coef_tr;
                crug = malha.coef_rug;
                min_t_color = malha.color;
            }
        }
    }

    for (int i = 0; i < n_spheres; i++){
        t = hit_sphere(s_list[i], r);
        if (t >= 1){
            if (t < min_t || min_t == -1 ){
                min_t = t;
                min_normal = sphere_normal(s_list[i].cent, t*r.direction() + r.origin());
                cd = s_list[i].coef_dif;
                ce = s_list[i].coef_esp;
                ca = s_list[i].coef_amb;
                cr = s_list[i].coef_ref;
                ct = s_list[i].coef_tr;
                crug = s_list[i].coef_rug;
                min_t_color = s_list[i].color;
            }
        }
    }

    for (int i = 0; i < n_planes; i++){
        t = hit_plane(p_list[i].pp, p_list[i].pv, r);
        if (t >= 1){
            if (t < min_t || min_t == -1 ){
                min_t = t;
                min_normal = p_list[i].pv;
                cd = p_list[i].coef_dif;
                ce = p_list[i].coef_esp;
                ca = p_list[i].coef_amb;
                cr = p_list[i].coef_ref;
                ct = p_list[i].coef_tr;
                crug = p_list[i].coef_rug;
                min_t_color = p_list[i].color;
            }
        }
    }

    point3 intersec = r.at(min_t);

    if (min_t != -1){
        min_t_color = phong_equation(ca, filter, l_list.size(), l_list, intersec,
                       min_t_color, cd, min_normal, ce, camera, crug);

    }

    // limitando valor rgb
    for (int i = 0; i < 3; i++){
        if (min_t_color.e[i] >= 255){
            min_t_color.e[i] = 255;
        }
    }

    return min_t_color;

}

point3 vec_matrix_mult(double** matrix_ptrs, point3 vec){
    double vec_sub[3] = {0., 0., 0.};

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            vec_sub[j] += vec[i] * matrix_ptrs[j][i];
        }
    }

    vec = point3(vec_sub[0], vec_sub[1], vec_sub[2]);
    return vec;
}   

int main() {

    // Image
    double aspect_ratio ;
    std::cin >> aspect_ratio;

    int image_width = 400;

    // Calculate the image height, and ensure that it's at least 1.
    int image_height = static_cast<int>(image_width / aspect_ratio);
    image_height = (image_height < 1) ? 1 : image_height;

    // Camera
    int x, y, z;
    std::cin >> x >> y >> z;
    auto camera_center = point3(x, y, z);
    std::cin >> x >> y >> z;
    auto v_up = vec3(x, y, z);
    std::cin >> x >> y >> z;
    auto M = point3(x, y, z);   // viewport center
    auto focal_length = (M-camera_center).length();   // distance between the viewport and the camera
    auto viewport_height = 2.0;
    auto viewport_width = viewport_height * (static_cast<double>(image_width)/image_height);

    // orthonormal vectors

    auto v_x = unit_vector(camera_center - M);
    v_x = unit_vector(v_x);
    auto v_z = -cross(v_up, -v_x);
    v_z = unit_vector(v_z);
    auto v_y = cross(v_z, v_x);
    v_y = unit_vector(v_y);
    

    // Calculate the vectors across the horizontal and down the vertical viewport edges.
    auto viewport_u = vec3(0, 0, viewport_width);
    auto viewport_v = vec3(0, -viewport_height, 0);  //top to bottom vector

    // Calculate the horizontal and vertical delta vectors from pixel to pixel.
    // little steps that go from a pixel's center to another pixel's center
    auto pixel_delta_u = viewport_u / image_width;
    auto pixel_delta_v = viewport_v / image_height;

    // Calculate the location of the upper left pixel.
    auto viewport_upper_left = M - viewport_u/2 - viewport_v/2; // top left corner
    auto top_left_pixel = viewport_upper_left + 0.5 * (pixel_delta_v + pixel_delta_u);  //  pixel (0, 0)

    // inputs para a malha de triangulos:
    int n_triangles, n_vertices;

    double matrix[3][3] = {{1.0, 0.0, 0.0  },
                           {0.0, 0.0, 1.0  },
                           {0.0, -1.0, 0.0 }};

    double matrix_B[3] = {0, 0, 0};

    double *matrix_ptrs[3] = {matrix[0], matrix[1], matrix[2]};

    std::cin >> n_triangles >> n_vertices;

    point3 *points_list = new point3[n_vertices];

    for (int i = 0; i < n_vertices; i++){
        std::cin >> x >> y >> z; 
        points_list[i] = point3(x,y,z);
    }

    bool af_transf = true; // ! true => aplica a transformaçao afim aos pontos da malha de triangulos

    if (af_transf){
        for (int i = 0; i < n_vertices; i++){
            points_list[i] = vec_matrix_mult(matrix_ptrs, points_list[i]);
            // aplicando a rotação

            for (int j = 0; j < 3; j++){
                points_list[i][j] += matrix_B[j];
            }
            // aplicando a translaçao

        }
    }

    point3 *triangles_list = new point3[n_triangles];
    for (int i = 0; i < n_triangles; i++){
        std::cin >> x >> y >> z;
        triangles_list[i] = point3(x, y, z); // se trata dos indices dos vertices que compoem cada triangulo
    }

    vec3 *triangles_normals = new vec3[n_triangles];
    for (int i = 0; i < n_triangles; i++){
        std::cin >> x >> y >> z;
        triangles_normals[i] = vec3(x, y, z); 
    }

    vec3 *vertices_normals = new vec3[n_vertices];
    for (int i = 0; i < n_vertices; i++){
        std::cin >> x >> y >> z;
        vertices_normals[i] = vec3(x, y, z); 
    }

    // criando objeto da malha de triangulos
    triangles malha = triangles(color(0, 0, 0), n_triangles, points_list, triangles_list,
                                0, 0, 0, 0, 0, 0);

    
    // definindo pontos de luz:
    std::vector<light> l_list;

    light l1 = light(point3(2, 5, 0), color(250, 250, 250));

    l_list.push_back(l1);

    color filter = color(10, 200, 10);

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    vec3 ray_direction, pixel_center, pixel_color;
    
    for (int j = 0; j < image_height; ++j) {
        for (int i = 0; i < image_width; ++i) {
            pixel_center = top_left_pixel + (i * pixel_delta_u) + (j * pixel_delta_v); // jumping from pixel to pixel with the mini-vectors
            ray_direction = pixel_center - camera_center;  // t == 1, as pixel_center = camera_center + t*ray_direction

            ray r(camera_center, ray_direction);

            pixel_color = ray_color(r, malha, l_list, filter, camera_center);   // painting the pixel with the color of the object that the pixel intercepted
            
            write_color(std::cout, pixel_color);
        }
    }

}