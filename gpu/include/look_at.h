#ifndef LOOK_AT_H
#define LOOK_AT_H

#include <Eigen/Core>
#include <Eigen/Geometry>

inline void look_at(const Eigen::Vector3f &eye, 
    const Eigen::Vector3f &up, 
    const Eigen::Vector3f &at,
    Eigen::Matrix4f &view);

inline void look_at(const Eigen::Vector3f &eye, 
    const Eigen::Vector3f &up, 
    const Eigen::Vector3f &at,
    Eigen::Matrix4f &view)
{
      //std::cout << eye << std::endl;
      Eigen::Vector3f dir = eye-at;
      Eigen::Vector3f Z = -dir.normalized();
      Eigen::Vector3f t = up;
      Eigen::Vector3f X = t.cross(Z);
      X.normalize();
      Eigen::Vector3f Y = Z.cross(X);
      Y.normalize();
      view.row(0) = Eigen::RowVector4f(X(0),X(1),X(2),-X.dot(eye));
      view.row(1) = Eigen::RowVector4f(Y(0),Y(1),Y(2),-Y.dot(eye));
      view.row(2) = Eigen::RowVector4f(Z(0),Z(1),Z(2),-Z.dot(eye));
      view.row(3) = Eigen::RowVector4f(0.,0.,0.,1.);

    //   view.transpose();
    //   Eigen::Affine3f transform(Eigen::Translation3f(-eye.head(3)));
    //   Eigen::Matrix4f translate = transform.matrix();
    //   view = translate * view;   
}

/*
vec3 zaxis = normal(eye - target);    // The "forward" vector.
    vec3 xaxis = normal(cross(up, zaxis));// The "right" vector.
    vec3 yaxis = cross(zaxis, xaxis);     // The "up" vector.

    // Create a 4x4 orientation matrix from the right, up, and forward vectors
    // This is transposed which is equivalent to performing an inverse 
    // if the matrix is orthonormalized (in this case, it is).
    mat4 orientation = {
       vec4( xaxis.x, yaxis.x, zaxis.x, 0 ),
       vec4( xaxis.y, yaxis.y, zaxis.y, 0 ),
       vec4( xaxis.z, yaxis.z, zaxis.z, 0 ),
       vec4(   0,       0,       0,     1 )
    };
    
    // Create a 4x4 translation matrix.
    // The eye position is negated which is equivalent
    // to the inverse of the translation matrix. 
    // T(v)^-1 == T(-v)
    mat4 translation = {
        vec4(   1,      0,      0,   0 ),
        vec4(   0,      1,      0,   0 ), 
        vec4(   0,      0,      1,   0 ),
        vec4(-eye.x, -eye.y, -eye.z, 1 )
    };

    // Combine the orientation and translation to compute 
    // the final view matrix. Note that the order of 
    // multiplication is reversed because the matrices
    // are already inverted.
    return ( orientation * translation );
*/

#endif