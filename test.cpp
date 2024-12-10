#include <stdio.h>
#include <string>

struct pp {
  float x;
  float y;
  float z;
  float w;
};

int main() {
    struct pp *ppos = new struct pp[10];
    ppos[0] = {1,2,3,4};
    ppos[1] = {10,20,30,40};
    ppos[2] = {11,21,31,41};
    ppos[3] = {12,22,32,42};

    float *x = (&ppos->z) + 4;
    printf("x,y,z,w = %f %f %f %f\n", ppos[0].x, ppos[0].y, ppos[0].z, ppos[0].w);
    printf("x %f\n", x[0]);
    // printf("x,y,z,w = %f %f %f %f\n", ppos[0].x, ppos[0].y, ppos[0].z, ppos[0].w);

    return 0;
}