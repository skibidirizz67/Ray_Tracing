#include <stdio.h>
#include <math.h>

typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    double r, g, b;
} Col;

typedef struct {
    Vec3 origin, direction;
} Ray;

Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 scalar_mult_vec3(double s, Vec3 v) {
    return (Vec3){v.x * s, v.y * s, v.z * s};
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

Vec3 at(Ray ray, double t) {
    return vec3_add(ray.origin, scalar_mult_vec3(t, ray.direction));
}

void print_progress(int bar_width, int progress, int total) {
    float ratio = (float)progress / total;
    int pos = bar_width * ratio;

    fprintf(stderr, "\r[");
    for (int i = 0; i < bar_width; i++) {
        if (i <= pos) {
            fprintf(stderr, "=");
        } else if (i == pos) {
            fprintf(stderr, ">");
        } else {
            fprintf(stderr, " ");
        }
    }
    fprintf(stderr, "] %d%% (%d/%d)", (int)(ratio * 100), progress, total);
}

void write_pixel(Vec3 pixel) {
    double r = pixel.x;
    double g = pixel.y;
    double b = pixel.z;

    int rb = (int)(255.999 * r);
    int gb = (int)(255.999 * g);
    int bb = (int)(255.999 * b);

    printf("%d %d %d\n", rb, gb, bb);
}

static inline Vec3 ray_color(Ray r) {
    Vec3 unit_direction = vec3_normalize(r.direction);
    double a = 0.5 * (unit_direction.y + 1.0);
    return vec3_add(
        scalar_mult_vec3((1.0 - a), (Vec3){1.0, 1.0, 1.0}),
        scalar_mult_vec3(a, (Vec3){0.5, 0.7, 1.0})
    );
}

int main() {
    double aspect_ratio  = 16.0 / 9.0;
    int img_width = 400;

    int img_height = (int)(img_width / aspect_ratio);
    if (img_height < 1) img_height = 1;

    double focal_len = 1.0;
    double viewport_height = 2.0;
    double viewport_width = viewport_height * ((double)img_width / img_height);
    Vec3 camera_center = (Vec3){0.0, 0.0, 0.0};

    Vec3 viewport_u = (Vec3){viewport_width, 0.0, 0.0};
    Vec3 viewport_v = (Vec3){0.0, -viewport_height, 0.0};

    Vec3 pixel_delta_u = vec3_div_scalar(viewport_u, img_width);
    Vec3 pixel_delta_v = vec3_div_scalar(viewport_v, img_height);

    Vec3 viewport_upper_left = vec3_sub(
        vec3_sub(
            vec3_sub(camera_center, (Vec3){0.0, 0.0, focal_len}),
            vec3_div_scalar(viewport_u, 2)
        ),
        vec3_div_scalar(viewport_v, 2)
    );
    Vec3 pixel00_loc = vec3_add(
        viewport_upper_left,
        scalar_mult_vec3(0.5, vec3_add(pixel_delta_u, pixel_delta_v))
    );

    printf("P3\n%d %d\n255\n", img_width, img_height);

    for (int y = 0; y < img_height; y++) {
        print_progress(20, y, img_height);
        for (int x = 0; x < img_width; x++) {
            Vec3 pixel_center = vec3_add(
                vec3_add(
                    pixel00_loc,
                    scalar_mult_vec3(x, pixel_delta_u)
                ),
                scalar_mult_vec3(y, pixel_delta_v)
            );
            Vec3 ray_direction = vec3_sub(pixel_center, camera_center);
            Ray r = (Ray){camera_center, ray_direction};

            Vec3 pixel_color = ray_color(r);
            write_pixel(pixel_color);
        }
    }

    fprintf(stderr, "\nDone.\n");
    
    return 0;
}