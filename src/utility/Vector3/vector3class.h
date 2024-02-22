#ifndef VECTOR3CLASS_H
#define VECTOR3CLASS_H

#include <vector>

template <typename T>
class Vector3 {
    
private:


public:
    T x,y,z;

    // Vector3();
    Vector3(T x0, T y0, T z0) {
        x=x0;
        y=y0;
        z=z0;
    };

    //void AddByVector(float,Vector3&);
    T norm();
    T d_square();
    T angle_with(Vector3&);
    T dot_with(Vector3&);


    // Vector3(std::vector<T>&);
    // Vector3(Vector3&,T);
    // Vector3(const Vector3&,const Vector3&);
    // Vector3(const Vector3&,const Vector3&,bool);

    inline T getX() const {return x;}
    inline T getY() const {return y;}
    inline T getZ() const {return z;}


    inline Vector3& operator=(const T &v2) {
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

    inline void operator+=(const T &f) {
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

    inline void operator-=(const T &f) {
        x -= f;
        y -= f;
        z -= f;
    }

    // v1 *= const
    inline void operator*=(const T &v2) {
        x *= v2;
        y *= v2;
        z *= v2;
    }

    // v1 /= const
    inline void operator/=(const T &v2) {
        x *= v2;
        y *= v2;
        z *= v2;
    }


};


#endif
