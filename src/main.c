#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define PI 3.1415926535897932385

#define SPHERE2H(obj) ((Hittable){(void*)&obj, &hit_sphere})

#define ASPECT_RATIO (16.0 / 9.0)
#define IMG_WIDTH 800


typedef struct {
    double x, y, z;
} Vec3;

Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 vec3_add_scalar(Vec3 v, double s) {
    return (Vec3){v.x + s, v.y + s, v.z + s};
}

Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 scalar_mult_vec3(double s, Vec3 v) {
    return (Vec3){v.x * s, v.y * s, v.z * s};
}

double vec3_len_pow2(Vec3 v) {
    return v.x*v.x + v.y*v.y + v.z*v.z;
}

double vec3_len(Vec3 v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vec3 vec3_normalize(Vec3 v) {
    Vec3 r = {0.0, 0.0, 0.0};
    double len = vec3_len(v);
    if (len > 0.0) {
        r.x = v.x / len;
        r.y = v.y / len;
        r.z = v.z / len;
    }
    return r;
}

Vec3 vec3_div_scalar(Vec3 v, double s) {
    Vec3 r = {0.0, 0.0, 0.0};
    if (s != 0.0) {
        r.x = v.x / s;
        r.y = v.y / s;
        r.z = v.z / s;
    }
    return r;
}

double vec3_dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}



typedef struct {
    double r, g, b;
} Col;

Col vec3_to_col(Vec3 v) {
    return (Col){v.x, v.y, v.z};
}



typedef struct {
    Vec3 center;
    double radius;
} Sphere;

typedef struct {
    Vec3 origin, direction;
} Ray;



typedef struct {
    double min, max;
} Range;

#define range_size(r) (r.max - r.min)
#define range_contains(r, x) (r.min <= x && x <= r.max)
#define range_surrounds(r, x) (r.min < x && x < r.max)



typedef struct {
    Vec3 p;
    Vec3 normal;
    double t;
    bool front_face;
} HitRecord;

typedef struct {
    void *data;
    bool (*hit)(void *data, Ray r, Range t, HitRecord *hrec);
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
    free(hlist);
}
bool hlist_hitall(HittableList *hlist, Ray r, Range t, HitRecord *hrec) {
    HitRecord temp_rec = {0};
    bool hit_any = false;
    double closest = t.max;

    for (size_t i = 0; i < hlist->size; i++) {
        if (hlist->hittables[i].data == NULL) continue;
        if (hlist->hittables[i].hit(hlist->hittables[i].data, r, (Range){t.min, closest}, &temp_rec)) {
            hrec = &temp_rec;
            hit_any = true;
            closest = temp_rec.t;
        }
    }

    return hit_any;
}



typedef struct {
    double aspect_ratio;
    int img_width, img_height;
    Vec3 center, pixel00_loc;
    Vec3 pixel_delta_u, pixel_delta_v;
} Camera;

Camera cam_init() {
    Camera cam = {0};

    cam.aspect_ratio = ASPECT_RATIO;
    cam.img_width = IMG_WIDTH;

    cam.img_height = (int)(cam.img_width / cam.aspect_ratio);
    if (cam.img_height < 1) cam.img_height = 1;

    cam.center = (Vec3){0.0, 0.0, 0.0};

    double focal_len = 1.0;
    double viewport_height = 2.0;
    double viewport_width = viewport_height * ((double)cam.img_width / cam.img_height);

    Vec3 viewport_u = (Vec3){viewport_width, 0.0, 0.0};
    Vec3 viewport_v = (Vec3){0.0, -viewport_height, 0.0};

    cam.pixel_delta_u = vec3_div_scalar(viewport_u, cam.img_width);
    cam.pixel_delta_v = vec3_div_scalar(viewport_v, cam.img_height);

    Vec3 viewport_upper_left = vec3_sub(
        vec3_sub(
            vec3_sub(cam.center, (Vec3){0.0, 0.0, focal_len}),
            vec3_div_scalar(viewport_u, 2)
        ),
        vec3_div_scalar(viewport_v, 2)
    );
    cam.pixel00_loc = vec3_add(
        viewport_upper_left,
        scalar_mult_vec3(0.5, vec3_add(cam.pixel_delta_u, cam.pixel_delta_v))
    );

    return cam;
}

Col ray_color(Ray r, HittableList *hlist) {
    HitRecord hrec;
    if (hlist_hitall(hlist, r, (Range){0, INFINITY}, &hrec)) {
        return vec3_to_col(scalar_mult_vec3(
            0.5,
            vec3_add_scalar(hrec.normal, 1)
        )); 
    }

    return (Col){0};
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

void write_pixel(Col pixel) {
    int rb = (int)(255.999 * pixel.r);
    int gb = (int)(255.999 * pixel.g);
    int bb = (int)(255.999 * pixel.b);

    printf("%d %d %d\n", rb, gb, bb);
}

void render(Camera *cam, HittableList *hlist) {
    printf("P3\n%d %d\n255\n", cam->img_width, cam->img_height);

    for (int y = 0; y < cam->img_height; y++) {
        print_progress(25, y, cam->img_height);
        for (int x = 0; x < cam->img_width; x++) {
            Vec3 pixel_center = vec3_add(
                vec3_add(
                    cam->pixel00_loc,
                    scalar_mult_vec3(x, cam->pixel_delta_u)
                ),
                scalar_mult_vec3(y, cam->pixel_delta_v)
            );
            Vec3 ray_direction = vec3_sub(pixel_center, cam->center);
            Ray r = (Ray){cam->center, ray_direction};

            Col pixel_color = ray_color(r, hlist);
            write_pixel(pixel_color);
        }
    }

    fprintf(stderr, "\nDone.\n");
}



double d2r(double deg) {
    return deg * PI / 180.0;
}

Vec3 at(Ray r, double t) {
    return vec3_add(r.origin, scalar_mult_vec3(t, r.direction));
}

void set_face_normal(HitRecord *hrec, Ray r, Vec3 outward_normal) {
    hrec->front_face = vec3_dot(r.direction, outward_normal) < 0;
    hrec->normal = hrec->front_face? outward_normal : scalar_mult_vec3(-1, outward_normal);
}

bool hit_sphere(void *data, Ray r, Range t, HitRecord *hrec) {
    Sphere *sphere = (Sphere*)data;

    double a = vec3_len_pow2(r.direction);
    Vec3 oc = vec3_sub(sphere->center, r.origin);
    double h = vec3_dot(r.direction, oc);
    double c = vec3_len_pow2(oc) - sphere->radius*sphere->radius;
    double discriminant = h*h - a*c;

    if (discriminant < 0) return false;
    
    double sqrtd = sqrt(discriminant);

    double root = (h - sqrtd) / a;
    if (!range_surrounds(t, root)) {
        root = (h + sqrtd) / a;
        if (!range_surrounds(t, root)) return false;
    }

    hrec->t = root;
    hrec->p = at(r, hrec->t);
    Vec3 outward_normal = vec3_div_scalar(vec3_sub(hrec->p, sphere->center), sphere->radius);
    set_face_normal(hrec, r, outward_normal);

    return true;
}

int main() {
    HittableList hlist = hlist_init();
    hlist_push(&hlist, SPHERE2H(((Sphere){(Vec3){0, 0, -1}, 0.5})));

    Camera cam = cam_init();
    render(&cam, &hlist);
    
    return 0;
}