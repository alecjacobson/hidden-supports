#ifndef PERSPECTIVE_PROJECTION_H
#define PERSPECTIVE_PROJECTION_H
#include "gl.h"
#include <Eigen/Core>
#include <iostream>

inline void perspective(float left, 
    float right, 
    float bottom, 
    float top, 
    float z_slice, 
    float far, 
    Eigen::Matrix4f &proj);

inline void perspective(float left, 
    float right, 
    float bottom, 
    float top, 
    float z_slice, 
    float far, 
    Eigen::Matrix4f &proj)
{
    proj(0,0) = 2.0*z_slice/(right-left);
    proj(1,1) = 2.0*z_slice/(top-bottom);
    proj(2,2) = -(far + z_slice)/(far - z_slice);
    proj(2,3) = -2.0*far*z_slice/(far - z_slice);
    proj(0,2) = (right+left)/(right-left);
    proj(3,2) = -1.0;
}

#endif