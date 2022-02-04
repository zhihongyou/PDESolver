//
// Created by xiuqi on 12/8/16.
//

#include <tgmath.h>
#include "Vector3.h"

Vector3::Vector3() : x(0.0),y(0.0),z(0.0){ }

Vector3::Vector3(float _x,float _y,float _z) : x(_x),y(_y),z(_z){ }

Vector3::Vector3(std::vector<float> &init) : x(init[0]),y(init[1]),z(init[2]){ }

Vector3::Vector3(Vector3 & init, float factor) : x(init.getX() * factor),y(init.getY() * factor),z(init.getZ() * factor){ }

Vector3::Vector3(const Vector3 & v1, const Vector3 & v2) : x(v1.getX() - v2.getX()), y(v1.getY() - v2.getY()), z(v1.getZ() - v2.getZ()){ }

Vector3::Vector3(const Vector3 &u, const Vector3 &v, bool): x(u.y*v.z - u.z*v.y),y(u.z*v.x - u.x * v.z),z(u.x* v.y - u.y*v.z) { }

/*
void Vector3::AddByVector(float scalar,Vector3 &v) {
    x += scalar * v.getX();
    y += scalar * v.getY();
    z += scalar * v.getZ();
}*/

float Vector3::norm() {
    return sqrt(d_square());
}

float Vector3::d_square() {
    return x*x + y*y + z*z;
}

float Vector3::angle_with(Vector3 &v) {
    float n_x,n_y,n_z;
    n_x = y * v.z - z * v.y;
    n_y = z * v.x - x * v.z;
    n_z = x * v.y - y * v.x;
    return atan(sqrt(n_x*n_x + n_y*n_y + n_z*n_z) / dot_with(v));
}

float Vector3::dot_with(Vector3 & v) {
    return x * v.getX() + y * v.getY() + z * v.getZ();
}

