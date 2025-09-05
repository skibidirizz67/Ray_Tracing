#include <stdio.h>

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

int main() {
    int img_width = 256;
    int img_height = 256;

    printf("P3\n%d %d\n255\n", img_width, img_height);

    for (int y = 0; y < img_height; y++) {
        print_progress(20, y, img_height);
        for (int x = 0; x < img_width; x++) {
            double r = (double)x / (img_width-1);
            double g = (double)y / (img_height-1);
            double b = 0.0;

            int ir = (int)(255.999 * r);
            int ig = (int)(255.999 * g);
            int ib = (int)(255.999 * b);

            printf("%d %d %d\n", ir, ig, ib);
        }
    }

    fprintf(stderr, "\nDone\n");
    
    return 0;
}