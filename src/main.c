#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926535897932385

#define SPHERE2H(obj) ((Hittable){(void*)&obj, &hit_sphere})
#define QUAD2H(obj) ((Hittable){(void*)&obj, &hit_quad})

#define LAMBERTIAN2M(mat) ((Material){(void*)&(mat), &scatter_lambertian, &emitted_default})
#define LAMBERTIANL2M(mat) ((Material){(void*)&(mat), &scatter_lambertian, &emitted_lambertian})
#define METAL2M(mat) ((Material){(void*)&(mat), &scatter_metal, &emitted_default})
#define DIELECTRIC2M(mat) ((Material){(void*)&(mat), &scatter_dielectric, &emitted_default})

static inline double dranddoub() { return rand() / (RAND_MAX + 1.0); }
static inline double randdoub(double min, double max) { return min + (max-min)*dranddoub(); }
static inline double d2r(double deg) { return deg * PI / 180.0; }
typedef struct {
    union {
        struct { double x, y, z; };
        struct { double r, g, b; };
    };
} Vec3;
static inline Vec3 vec3_add(Vec3 a, Vec3 b) { return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline Vec3 vec3_add_scalar(Vec3 v, double s) { return (Vec3){v.x + s, v.y + s, v.z + s}; }
static inline Vec3 vec3_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline Vec3 vec3_hadamard(Vec3 a, Vec3 b) { return (Vec3){a.x*b.x, a.y*b.y, a.z*b.z}; }
static inline double vec3_dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 vec3_cross(Vec3 a, Vec3 b) { return (Vec3){a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
static inline Vec3 scalar_mult_vec3(double s, Vec3 v) { return (Vec3){v.x * s, v.y * s, v.z * s}; }
static inline Vec3 vec3_div_scalar(Vec3 v, double s) { return (s==0.0)? (Vec3){0} : (Vec3){v.x/s, v.y/s, v.z/s}; }
static inline double vec3_len_pow2(Vec3 v) { return v.x*v.x + v.y*v.y + v.z*v.z; }
static inline double vec3_len(Vec3 v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
static inline Vec3 vec3_normalize(Vec3 v) {
    Vec3 r = {0.0, 0.0, 0.0};
    double len = vec3_len(v);
    if (len > 0.0) {
        r.x = v.x / len;
        r.y = v.y / len;
        r.z = v.z / len;
    }
    return r;
}
static inline bool vec3_near_zero(Vec3 v) { return (fabs(v.x) < 1e-8) && (fabs(v.y) < 1e-8) && (fabs(v.z) < 1e-8); }
static inline Vec3 vec3_random(double min, double max) { return (Vec3){randdoub(min,max), randdoub(min,max), randdoub(min,max)}; }
static inline Vec3 vec3_drandom() { return (Vec3){dranddoub(), dranddoub(), dranddoub()}; }
static inline Vec3 vec3_rand_norm() {
    for (;;) {
        Vec3 p = { randdoub(-1,1), randdoub(-1,1), randdoub(-1,1) };
        double len2 = vec3_len_pow2(p);
        if (len2 > 1e-12 && len2 <= 1.0) return vec3_div_scalar(p, sqrt(len2));
    }
}
static inline Vec3 vec3_random_on_hemisphere(Vec3 normal) {
    Vec3 v = vec3_rand_norm();
    return (vec3_dot(v, normal) > 0.0)? v : scalar_mult_vec3(-1.0, v);
}
static inline Vec3 vec3_reflect(Vec3 v, Vec3 n) { return vec3_sub(v, scalar_mult_vec3(2*vec3_dot(v, n), n)); }
static inline Vec3 vec3_refract(Vec3 uv, Vec3 n, double ri) {
    double cos_theta = fmin(vec3_dot(scalar_mult_vec3(-1.0, uv), n), 1.0);
    Vec3 r_perp = scalar_mult_vec3(ri, vec3_add(uv, scalar_mult_vec3(cos_theta, n)));
    double len2 = vec3_len_pow2(r_perp);
    double under = 1.0 - len2;
    double root = (under > 0.0) ? -sqrt(under) : 0.0;
    Vec3 r_par = scalar_mult_vec3(root, n);
    return vec3_add(r_perp, r_par);
}

static inline double clamp(double n, double min, double max) { return (n < min)? min : ((n > max)? max : n); }



typedef struct {
    Vec3 p;
    Vec3 normal;
    double t;
    bool front_face;
} HitRecord;

typedef struct {
    Vec3 origin, direction;
} Ray;

static inline Vec3 at(Ray r, double t) {
    return vec3_add(r.origin, scalar_mult_vec3(t, r.direction));
}

typedef struct {
	void *data;
	bool (*scatter)(void *data, Ray r, HitRecord *hrec, Vec3 *atten, Ray *scattered);
    Vec3 (*emitted)(void *data);
} Material;

typedef struct {
    void *data;
    bool (*hit)(void *data, Ray r, double tmin, double tmax, HitRecord *hrec, Material *mat);
} Hittable;
typedef struct {
    Hittable *hittables;
    size_t size, capacity;
} HittableList;
HittableList hlist_init() {
    HittableList hlist = {0};
    hlist.hittables = (Hittable*)calloc(1, sizeof(Hittable));
    if (hlist.hittables == NULL) {
        fprintf(stderr, "Failed to init HittableList");
        return (HittableList){0};
    }
    hlist.capacity = 1;
    return hlist;
}
void hlist_resize(HittableList *hlist) {
    Hittable *temp = realloc(hlist->hittables, sizeof(Hittable)*(hlist->capacity*2));
    if (temp == NULL) {
        fprintf(stderr, "Failed to resize HittableList");
        return;
    }
    hlist->hittables = temp;
    hlist->capacity *= 2;
}
void hlist_push(HittableList *hlist, Hittable hittable) {
    if (hlist->size == hlist->capacity) hlist_resize(hlist);
    hlist->hittables[hlist->size] = hittable;
    hlist->size++;
}
void hlist_pop(HittableList *hlist) {
    hlist->hittables[hlist->size-1] = (Hittable){0};
    hlist->size--;
}
void hlist_clear(HittableList *hlist) {
    for (size_t i = 0; i < hlist->size; i++) hlist->hittables[i] = (Hittable){0};
    hlist->size = 0;
}
void hlist_destroy(HittableList *hlist) {
    free(hlist->hittables);
    hlist = (HittableList*){0};
}
bool hlist_hitall(HittableList *hlist, Ray r, double tmin, double tmax, HitRecord *hrec, Material *mat) {
    HitRecord temp_rec = {0};
    bool hit_any = false;
    double closest = tmax;

    for (size_t i = 0; i < hlist->size; i++) {
        if (hlist->hittables[i].data == NULL) continue;
        if (hlist->hittables[i].hit(hlist->hittables[i].data, r, tmin, closest, &temp_rec, mat)) {
            *hrec = temp_rec;
            hit_any = true;
            closest = temp_rec.t;
        }
    }

    return hit_any;
}


typedef struct {
	Vec3 albedo;
} Lambertian;
bool scatter_lambertian(void *data, Ray r, HitRecord *hrec, Vec3 *attenuation, Ray *scattered) {
	Lambertian *this = (Lambertian*)data;
	Vec3 dir = vec3_add(hrec->normal, vec3_rand_norm());
	if (vec3_near_zero(dir)) dir = hrec->normal;
	*scattered = (Ray){hrec->p, dir};
	*attenuation = this->albedo;
	return true; 
}

typedef struct {
	Vec3 albedo;
    double roughness;
} Metal;
bool scatter_metal(void *data, Ray r, HitRecord *hrec, Vec3 *attenuation, Ray *scattered) {
	Metal *this = (Metal*)data;
	Vec3 reflected = vec3_reflect(r.direction, hrec->normal);
    reflected = vec3_add(vec3_normalize(reflected), scalar_mult_vec3(this->roughness, vec3_rand_norm()));
	*scattered = (Ray){hrec->p, reflected};
	*attenuation = this->albedo;
	return vec3_dot(scattered->direction, hrec->normal) > 0;
}

typedef struct {
	double refraction_index;
} Dielectric;
double rSchlick(Vec3 incident, Vec3 normal, double n1, double n2) {
    double r0 = (n1 - n2) / (n1 + n2);
    r0 *= r0;
    
    float x = 1.0 + vec3_dot(normal, incident);
    return clamp(r0 + (1.0-r0)*x*x*x*x*x, 0, 1);
}
// #TODO: blend two rays instead of picking random
bool scatter_dielectric(void *data, Ray r, HitRecord *hrec, Vec3 *attenuation, Ray *scattered) {
    Dielectric *this = (Dielectric*)data;
    *attenuation = (Vec3){1.0, 1.0, 1.0};
    double n1 = hrec->front_face? 1.0 : this->refraction_index;
    double n2 = hrec->front_face? this->refraction_index : 1.0;
    double eta = n1 / n2;
    Vec3 unit_direction = vec3_normalize(r.direction);
    double cos_theta = fmin(vec3_dot(
        scalar_mult_vec3(-1, unit_direction),
        hrec->normal
    ), 1.0);
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    bool cant_refract = (eta * sin_theta) > 1.0;
    Vec3 direction;
    double reflect_prob = rSchlick(unit_direction, hrec->normal, n1, n2);
    if (cant_refract || reflect_prob > dranddoub()) direction = vec3_reflect(unit_direction, hrec->normal);
    else direction = vec3_refract(unit_direction, hrec->normal, eta);
    *scattered = (Ray){hrec->p, direction};
    return true;
}


typedef struct {
    Vec3 albedo;
} LambertianLight;
Vec3 emitted_default(void *data) {
    return (Vec3){0.0, 0.0, 0.0};
}
Vec3 emitted_lambertian(void *data) {
    LambertianLight *this = (LambertianLight*)data;
    return this->albedo;
}


typedef struct {
    Vec3 center;
    double radius;
    Material mat;
} Sphere;
static inline void set_face_normal(HitRecord *hrec, Ray r, Vec3 outward_normal) {
    hrec->front_face = vec3_dot(r.direction, outward_normal) < 0;
    hrec->normal = hrec->front_face? outward_normal : scalar_mult_vec3(-1, outward_normal);
}
bool hit_sphere(void *data, Ray r, double tmin, double tmax, HitRecord *hrec, Material *mat) {
    Sphere *this = (Sphere*)data;

    double a = vec3_len_pow2(r.direction);
    Vec3 oc = vec3_sub(this->center, r.origin);
    double h = vec3_dot(r.direction, oc);
    double c = vec3_len_pow2(oc) - this->radius*this->radius;
    double discriminant = h*h - a*c;

    if (discriminant < 0) return false;
    
    double sqrtd = sqrt(discriminant);

    double root = (h - sqrtd) / a;
    if (!(tmin < root && root < tmax)) {
        root = (h + sqrtd) / a;
        if (!(tmin < root && root < tmax)) return false;
    }

    hrec->t = root;
    hrec->p = at(r, hrec->t);
    Vec3 outward_normal = vec3_div_scalar(vec3_sub(hrec->p, this->center), this->radius);
    set_face_normal(hrec, r, outward_normal);
    *mat = this->mat;

    return true;
}

typedef struct {
    Vec3 Q;
    Vec3 u, v, w;
    Vec3 normal;
    double D;
    Material mat;
} Quad;
void quad_init(Quad *quad) {
    Vec3 n = vec3_cross(quad->u, quad->v);
    quad->normal = vec3_normalize(n);
    quad->D = vec3_dot(quad->normal, quad->Q);
    quad->w = scalar_mult_vec3(1.0/vec3_dot(n, n), n);
}
bool is_interior(double a, double b, HitRecord *hrec) {
    if (!(a >= 0 && a <= 1) || !(b >= 0 && b <= 1)) return false;

    //hrec->u = a;
    //hrec->v = b;
    return true;
}
bool hit_quad(void *data, Ray r, double tmin, double tmax, HitRecord *hrec, Material *mat) {
    Quad *this = (Quad*)data;

    double denom = vec3_dot(this->normal, r.direction);

    if (fabs(denom) < 1e-8) return false;

    double t = (this->D - vec3_dot(this->normal, r.origin)) / denom;
    if (!(t >= tmin && t <= tmax)) return false;

    Vec3 intersection = at(r, t);
    Vec3 planar_hitpt_vector = vec3_sub(intersection, this->Q);
    double alpha = vec3_dot(this->w, vec3_cross(planar_hitpt_vector, this->v));
    double beta = vec3_dot(this->w, vec3_cross(this->u, planar_hitpt_vector));

    if (!is_interior(alpha, beta, hrec)) return false;

    hrec->t = t;
    hrec->p = intersection;
    *mat = this->mat;
    set_face_normal(hrec, r, this->normal);

    return true;
}


typedef struct {
    double aspect_ratio, vfov;
    int img_width, img_height;
    int samples_per_pixel, max_depth;
    double pixel_samples_scale;
    Vec3 lookfrom, lookat, lookup;
    Vec3 u, v, w;
    Vec3 pixel00_loc;
    Vec3 pixel_delta_u, pixel_delta_v;
    Vec3 background;
} Camera;

void cam_init(Camera *cam) {
    cam->img_height = (int)(cam->img_width / cam->aspect_ratio);
    if (cam->img_height < 1) cam->img_height = 1;

    cam->pixel_samples_scale = 1.0 / cam->samples_per_pixel;

    double focal_len = vec3_len(vec3_sub(cam->lookfrom, cam->lookat));
    double viewport_height = 2.0 * tan(d2r(cam->vfov)/2.0) * focal_len;
    double viewport_width = viewport_height * ((double)cam->img_width / cam->img_height);

    cam->w = vec3_normalize(vec3_sub(cam->lookfrom, cam->lookat));
    cam->u = vec3_normalize(vec3_cross(cam->lookup, cam->w));
    cam->v = vec3_cross(cam->w, cam->u);

    Vec3 viewport_u = scalar_mult_vec3(viewport_width, cam->u);
    Vec3 viewport_v = scalar_mult_vec3(viewport_height, cam->v);

    cam->pixel_delta_u = vec3_div_scalar(viewport_u, cam->img_width);
    cam->pixel_delta_v = vec3_div_scalar(viewport_v, cam->img_height);

    Vec3 viewport_upper_left = vec3_sub(
        vec3_sub(
            vec3_sub(cam->lookfrom, scalar_mult_vec3(focal_len, cam->w)),
            vec3_div_scalar(viewport_u, 2)
        ),
        vec3_div_scalar(viewport_v, 2)
    );
    cam->pixel00_loc = vec3_add(
        viewport_upper_left,
        scalar_mult_vec3(0.5, vec3_add(cam->pixel_delta_u, cam->pixel_delta_v))
    );
}

Vec3 ray_color(Ray r, int depth, HittableList *hlist, Vec3 background) {
    if (depth <= 0) return (Vec3){0.0, 0.0, 0.0};

    HitRecord hrec = {0};
    Material mat;
    if (!hlist_hitall(hlist, r, 0.001, INFINITY, &hrec, &mat)) return background;

    Ray scattered;
    Vec3 attenuation;
    Vec3 color_emission = mat.emitted(mat.data);
    if (!mat.scatter(mat.data, r, &hrec, &attenuation, &scattered)) return color_emission;

    Vec3 color_scatter = vec3_hadamard(attenuation, ray_color(scattered, depth-1, hlist, background));

    return vec3_add(color_emission, color_scatter);
}

void print_progress(int bar_width, int progress, int total) {
    double ratio = (double)progress / total;
    int pos = bar_width * ratio;

    fprintf(stderr, "\r[");
    for (int i = 0; i < bar_width; i++) {
        if (i < pos)
            fprintf(stderr, "=");
        else if (i == pos)
            fprintf(stderr, ">");
        else
            fprintf(stderr, " ");
    }
    fprintf(stderr, "] %d%% (%d/%d)", (int)(ratio * 100), progress, total);
}

double linear_to_gamma(double linear) {
    if (linear > 0) return sqrt(linear);
    return 0;
}
void write_pixel(Vec3 pixel) {   
    double r = linear_to_gamma(pixel.r);
    double g = linear_to_gamma(pixel.g);
    double b = linear_to_gamma(pixel.b);

    int rb = (int)(256 * clamp(r, 0, 0.999));
    int gb = (int)(256 * clamp(g, 0, 0.999));
    int bb = (int)(256 * clamp(b, 0, 0.999));

    printf("%d %d %d\n", rb, gb, bb);
}

Vec3 sample_square() {
    return (Vec3){randdoub(0, 1)-0.5, randdoub(0, 1)-0.5, 0.0};
}

Ray get_ray(Camera *cam, int x, int y) {
    Vec3 offset = sample_square();
    Vec3 pixel_sample = vec3_add(
        vec3_add(
            cam->pixel00_loc,
            scalar_mult_vec3(x + offset.x, cam->pixel_delta_u)
        ),
        scalar_mult_vec3(y + offset.y, cam->pixel_delta_v)
    );
    Vec3 direction = vec3_sub(pixel_sample, cam->lookfrom);
    return (Ray){cam->lookfrom, direction};
}

void render(Camera *cam, HittableList *hlist) {
    fprintf(stderr, "Rendering %dx%d | %d samples | %d depth\n", cam->img_width, cam->img_height, cam->samples_per_pixel, cam->max_depth);

    printf("P3\n%d %d\n255\n", cam->img_width, cam->img_height);
    for (int y = 0; y < cam->img_height; y++) {
        print_progress(25, y, cam->img_height);
        for (int x = 0; x < cam->img_width; x++) {
            Vec3 pixel_color = {0};
            for (int sample = 0; sample < cam->samples_per_pixel; sample++) {
                Ray r = get_ray(cam, x, y);
                pixel_color = vec3_add(pixel_color, ray_color(r, cam->max_depth, hlist, cam->background));
            }
            write_pixel(scalar_mult_vec3(cam->pixel_samples_scale, pixel_color));
        }
    }

    fprintf(stderr, "\nDone.\n");
}

void scene0() {
    Camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.vfov = 60;
    cam.img_width = 480;
    cam.samples_per_pixel = 32;
    cam.max_depth = 6;
    cam.lookfrom = (Vec3){0, 0, 0};
    cam.lookat = (Vec3){0, 0, -1};
    cam.lookup = (Vec3){0, -1, 0};
    cam.background = (Vec3){0.7, 0.8, 1.0};
    cam_init(&cam);

    Lambertian diffr = {(Vec3){1, 0, 0}};
	Lambertian diffb = {(Vec3){0, 0, 1}};
	Lambertian diffg = {(Vec3){0, 1, 0}};

    Dielectric glass = {1.33};
    Dielectric air = {1/1.33};

	Metal met = {(Vec3){0.8, 0.8, 0.8}, 0.1};

	Sphere obj0 = {
		(Vec3){-1, 0, -1}, 0.5, LAMBERTIAN2M(diffr)
	};
	Sphere obj1 = {
		(Vec3){0, -100.5, -1}, 100, LAMBERTIAN2M(diffb)
	};
    Sphere obj2 = {
		(Vec3){0, 0.25, -1.4}, 0.5, METAL2M(met)
	};
	Sphere obj3 = {
		(Vec3){1, 0, -1}, 0.5, DIELECTRIC2M(glass)
	};
    Sphere obj4 = {
		(Vec3){1, 0, -1}, 0.4, DIELECTRIC2M(air)
	};
    Quad obj6 = {
		(Vec3){-1, 0, -2}, (Vec3){3, 0, -2}, (Vec3){1, 4, -2}, .mat=LAMBERTIAN2M(diffg)
	};
    quad_init(&obj6);

    HittableList hlist = hlist_init();
    hlist_push(&hlist, SPHERE2H(obj0));
    hlist_push(&hlist, SPHERE2H(obj1));
    hlist_push(&hlist, SPHERE2H(obj2));
    hlist_push(&hlist, SPHERE2H(obj3));
    hlist_push(&hlist, SPHERE2H(obj4));
    hlist_push(&hlist, QUAD2H(obj6));

    render(&cam, &hlist);
    hlist_destroy(&hlist);
}

void scene_light() {
    Camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.vfov = 60;
    cam.img_width = 480;
    cam.samples_per_pixel = 64;
    cam.max_depth = 64;
    cam.lookfrom = (Vec3){0, 0, 0};
    cam.lookat = (Vec3){0, 0, -1};
    cam.lookup = (Vec3){0, -1, 0};
    cam.background = (Vec3){0, 0, 0};
    cam_init(&cam);

    Lambertian diffr = {(Vec3){1, 0, 0}};
	Lambertian diffb = {(Vec3){0, 0, 1}};

    Dielectric glass = {1.33};
    Dielectric air = {1/1.33};

	Metal met = {(Vec3){0.8, 0.8, 0.8}, 0.1};

    LambertianLight light = {(Vec3){15.0, 15.0, 15.0}};

	Sphere obj0 = {
		(Vec3){-1, 0, -1}, 0.5, LAMBERTIAN2M(diffr)
	};
	Sphere obj1 = {
		(Vec3){0, -100.5, -1}, 100, LAMBERTIAN2M(diffb)
	};
    Sphere obj2 = {
		(Vec3){0, 0.25, -1.4}, 0.5, METAL2M(met)
	};
	Sphere obj3 = {
		(Vec3){1, 0, -1}, 0.5, DIELECTRIC2M(glass)
	};
    Sphere obj4 = {
		(Vec3){1, 0, -1}, 0.4, DIELECTRIC2M(air)
	};
    Sphere obj5 = {
		(Vec3){0.5, 1.0, 0}, 0.25, LAMBERTIANL2M(light)
	};

    HittableList hlist = hlist_init();
    hlist_push(&hlist, SPHERE2H(obj0));
    hlist_push(&hlist, SPHERE2H(obj1));
    hlist_push(&hlist, SPHERE2H(obj2));
    hlist_push(&hlist, SPHERE2H(obj3));
    hlist_push(&hlist, SPHERE2H(obj4));
    hlist_push(&hlist, SPHERE2H(obj5));

    render(&cam, &hlist);
    
    hlist_destroy(&hlist);
}

void cornell_box() {
    Camera cam;
    cam.aspect_ratio = 1.0;
    cam.vfov = 39;
    cam.img_width = 180;
    cam.samples_per_pixel = 50;
    cam.max_depth = 64;
    cam.background = (Vec3){0, 0, 0};
    cam.lookfrom = (Vec3){278, 278, -800};
    cam.lookat = (Vec3){278, 278, 0};
    cam.lookup = (Vec3){0, -1, 0};
    cam_init(&cam);

    Material red = LAMBERTIAN2M(((Lambertian){(Vec3){0.65, 0.05, 0.05}}));
    Material white = LAMBERTIAN2M(((Lambertian){(Vec3){0.73, 0.73, 0.73}}));
    Material green = LAMBERTIAN2M(((Lambertian){(Vec3){0.12, 0.45, 0.15}}));
    Material light = LAMBERTIANL2M(((LambertianLight){(Vec3){15, 15, 15}}));

    Quad q0 = {(Vec3){555, 0, 0}, (Vec3){0, 555, 0}, (Vec3){0, 0, 555}, .mat=green};
    Quad q1 = {(Vec3){0, 0, 0}, (Vec3){0, 555, 0}, (Vec3){0, 0, 555}, .mat=red};
    Quad q2 = {(Vec3){343, 554, 332}, (Vec3){-130, 0, 0}, (Vec3){0, 0, -105}, .mat=light};
    Quad q3 = {(Vec3){0, 0, 0}, (Vec3){555, 0, 0}, (Vec3){0, 0, 555}, .mat=white};
    Quad q4 = {(Vec3){555, 555, 555}, (Vec3){-555, 0, 0}, (Vec3){0, 0, -555}, .mat=white};
    Quad q5 = {(Vec3){0, 0, 555}, (Vec3){555, 0, 0}, (Vec3){0, 555, 0}, .mat=white};
    quad_init(&q0);
    quad_init(&q1);
    quad_init(&q2);
    quad_init(&q3);
    quad_init(&q4);
    quad_init(&q5);

    HittableList hlist = hlist_init();
    hlist_push(&hlist, QUAD2H(q0));
    hlist_push(&hlist, QUAD2H(q1));
    hlist_push(&hlist, QUAD2H(q2));
    hlist_push(&hlist, QUAD2H(q3));
    hlist_push(&hlist, QUAD2H(q4));
    hlist_push(&hlist, QUAD2H(q5));

    render(&cam, &hlist);
    hlist_destroy(&hlist);
}

void glass_sphere() {
    Camera cam;
    cam.aspect_ratio = 16.0 / 9.0;
    cam.vfov = 40;
    cam.img_width = 600;
    cam.samples_per_pixel = 128;
    cam.max_depth = 6;
    cam.lookfrom = (Vec3){0, 0, 1};
    cam.lookat = (Vec3){0, 0, 0};
    cam.lookup = (Vec3){0, -1, 0};
    cam.background = (Vec3){0.0, 0.0, 0.0};
    cam_init(&cam);

    Lambertian ground = {(Vec3){0.0, 0.0, 0.9}};
    Lambertian red = {(Vec3){0.9, 0.0, 0.0}};
    Dielectric glass = {1.33};
    LambertianLight light = {(Vec3){15.0, 15.0, 15.0}};

    Sphere obj0 = {
		(Vec3){0, -100.5, -1}, 100, LAMBERTIAN2M(ground)
	};
    Sphere obj1 = {
        (Vec3){0, 0, -1}, 0.5, DIELECTRIC2M(glass)
    };
    Quad obj2 = {
        (Vec3){3, 2, -1}, (Vec3){0, 0, -1}, (Vec3){1, -1, 0}, .mat=LAMBERTIANL2M(light)
    };
    quad_init(&obj2);

    HittableList hlist = hlist_init();
    hlist_push(&hlist, SPHERE2H(obj0));
    hlist_push(&hlist, SPHERE2H(obj1));
    hlist_push(&hlist, QUAD2H(obj2));

    render(&cam, &hlist);
    hlist_destroy(&hlist);
}

int main() {
    srand(time(NULL));

    scene0();
    
    return 0;
}