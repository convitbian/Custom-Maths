#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <CustomMath.h>

struct Vector {
    float x;
    float y;
    
    Vector(float vec_x, float vec_y) {
        x = vec_x;
        y = vec_y;
    }

    Vector add(Vector add_vec) {
        return Vector(x + add_vec.x, y + add_vec.y);
    }

    Vector substract(Vector sub_vec) {
        return Vector(x - sub_vec.x, y - sub_vec.y);
    }

    float magnitude() {
        return c_sqrt(x * x + y * y);
    }

    Vector normal() {
        return Vector(x / magnitude(), y / magnitude());
    }

    Vector scale(float scaler) {
        return Vector(x * scaler, y * scaler);
    }

    Vector negate() {
        return Vector(-x, -y);
    }

    Vector perpendicular() {  // up to left
        return Vector(-y, x);
    }

    Vector rotate(float degree) {
        // x' = x * cos(a) - y * sin(a)
        // y' = x * sin(a) + y * cos(a)

        float x_new = x * c_cos(degree) - y * c_sin(degree);
        float y_new = x * c_sin(degree) + y * c_cos(degree);
        return Vector(x_new, y_new);
    }

    float dot(Vector vecA, Vector vecB) {
        return vecA.x * vecB.x + vecA.y * vecB.y;
    }

    float cross(Vector vecA, Vector vecB) {
        return vecA.x * vecB.y - vecA.y * vecB.x;
    }

    float deg_vectors(Vector vecA, Vector vecB) {
        float dot_value = dot(vecA, vecB);
        float mag_squared = vecA.magnitude() * vecB.magnitude();
        return c_arccos(dot_value / mag_squared);
    }

    float distance(Vector endpoint1, Vector endpoint2, Vector point, bool square_for_speed = false) {
        float a = (endpoint2.y - endpoint1.y) / (endpoint2.x - endpoint1.x);
        float c = endpoint1.y - a * endpoint1.x;

        float distance_squared = c_pow(a * point.x - point.y + c, 2) / (c_pow(a, 2) + 1);

        if(square_for_speed == true) return distance_squared;
        else return c_sqrt(distance_squared);

        return 0;
    }
};

#endif
