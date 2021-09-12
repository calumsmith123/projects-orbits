/*
 * Structure to store position and motion of stars, plantes, moons etc
 */

#define MAX_NAME 32u

typedef enum Dimensions {
    X, Y, Z, N
} dimension;

struct Body {
    char name[MAX_NAME];
    double mass;
    double radius[N];
    double velocity[N];
    double acceleration[N];
} body;
