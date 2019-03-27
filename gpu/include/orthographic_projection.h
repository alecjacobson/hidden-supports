#ifndef ORTHOGRAPHIC_PROJECTION_H
#define ORTHOGRAPHIC_PROJECTION_H
#include "gl.h"
#include <Eigen/Core>
#include <iostream>

inline void orthographic(float left, 
    float right, 
    float bottom, 
    float top, 
    float z_slice, 
    float far, 
    Eigen::Matrix4f &proj);

inline void orthographic(float left, 
    float right, 
    float bottom, 
    float top, 
    float z_slice, 
    float far, 
    Eigen::Matrix4f &proj)
{
    proj(0,0) = 2.0/(right - left);
    proj(1,1) = 2.0/(top - bottom);
    proj(2,2) = 2.0/(z_slice - far);
    proj(3,3) = 1.0;
    proj(0,3) = -(right + left)/(right - left);;
    proj(1,3) = -(top + bottom)/(top - bottom);
    proj(2,3) = -(far + z_slice)/(far - z_slice);
}

#endif