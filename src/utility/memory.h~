//
// Created by xiuqi on 12/8/16.
//

#ifndef READ_PDB_VECTOR3_H
#define READ_PDB_VECTOR3_H

#include <vector>


class Vector3 {
private:


public:
    float x,y,z;

    //void AddByVector(float,Vector3&);
    float norm();
    float d_square();
    float angle_with(Vector3&);
    float dot_with(Vector3&);


    Vector3();
    Vector3(float,float,float);
    Vector3(std::vector<float>&);
    Vector3(Vector3&,float);
    Vector3(const Vector3&,const Vector3&);

    Vector3(const Vector3&,const Vector3&,bool);

    inline float getX() const {return x;}
    inline float getY() const {return y;}
    inline float getZ() const {return z;}


    inline Vector3& operator=(const float &v2) {
        x = v2;
        y = v2;
        z = v2;
        return *this;
    }

    //  v1 += v2;
    inline void operator+=(const Vector3 &v2) {
        x += v2.x;
        y += v2.y;
        z += v2.z;
    }

    inline void operator+=(const float &f) {
        x += f;
        y += f;
        z += f;
    }

    // v1 -= v2;
    inline void operator-=(const Vector3 &v2) {
        x -= v2.x;
        y -= v2.y;
        z -= v2.z;
    }

    inline void operator-=(const float &f) {
        x -= f;
        y -= f;
        z -= f;
    }

    // v1 *= const
    inline void operator*=(const float &v2) {
        x *= v2;
        y *= v2;
        z *= v2;
    }

    // v1 /= const
    inline void operator/=(const float& v2) {
        x *= v2;
        y *= v2;
        z *= v2;
    }






};


#endif //READ_PDB_VECTOR3_H
