#include <iostream>
#include "vec3.h"
#include "ray.h"
#include "color.h"
#include "sphere.h"
#include "plane.h"

#define EPSILON 2.22045e-016  // menor diferença entre dois doubles

double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
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

color ray_color(const ray& r, point3* points_list, point3* triangles_list, int n_triangles) {

    int n_objects = 3;  // numero de objetos que serao renderizados, fora os triangulos
    double min_t = -1;
    color min_t_color = color(0, 0, 0);

    double *o_list = new double[n_objects];
    double *t_list = new double[n_triangles];

    color *color_list = new color[n_objects];

    
    o_list[0] = hit_sphere(point3(2,0,0), 0.2, r);
    o_list[1] = hit_sphere(point3(5,5,0), 0.4, r);
    o_list[2] = hit_plane(point3(7, 2, 2), vec3(1, 0.2, 0.3), r);

    for(int i = 0; i < n_triangles; i++){
        t_list[i] = hit_triangle(triangles_list[i], points_list, r);
    }

    color_list[0] = color(1, 1, 1);
    color_list[1] = color(1, 0, 0);
    color_list[2] = color(0, 0, 1);

    color t_color = color(1,1,0);


    for (int i = 0; i < n_objects; i++){
        if (o_list[i] >= 1){
            if (o_list[i] < min_t || min_t == -1 ){
                min_t = o_list[i];
                min_t_color = color_list[i];
            }
        }
    }

    for (int i = 0; i < n_triangles; i++){
        if (t_list[i] >= 1){
            if (t_list[i] < min_t || min_t == -1){
                min_t = t_list[i];
                min_t_color = t_color;
            }
        }
    }

    return min_t_color;

}

/*
    A = [ cos(a) sen(a) 0]
        [-sen(a) cos(a) 0]
        [ 0      0      1]

        matriz de rotaçao em torno do eixo z

    B = [0]
        [1]
        [1]

        translaçao


                 [x*cos(a) +  y*sen(a)]     [0]
    T(x, y, z) = [-x*sen(a) + y*cos(a)]  +  [1]
                 [          z         ]     [1]
*/

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

    

    // Render

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    vec3 ray_direction, pixel_center, pixel_color;
    
    for (int j = 0; j < image_height; ++j) {
        for (int i = 0; i < image_width; ++i) {
            pixel_center = top_left_pixel + (i * pixel_delta_u) + (j * pixel_delta_v); // jumping from pixel to pixel with the mini-vectors
            ray_direction = pixel_center - camera_center;  // t == 1, as pixel_center = camera_center + t*ray_direction

            ray r(camera_center, ray_direction);

            pixel_color = ray_color(r, points_list, triangles_list, n_triangles);   // painting the pixel with the color of the object that the pixel intercepted
            
            write_color(std::cout, pixel_color);
        }
    }

}