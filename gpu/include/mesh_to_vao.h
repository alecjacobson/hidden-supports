#ifndef MESH_TO_VAO_H
#define MESH_TO_VAO_H
#include "gl.h"
#include <Eigen/Core>

// Send a triangle mesh to the GPU using a vertex array object.
//
// Inputs:
//   V  #V by 3 list of 3D mesh vertex positions
//   F  #F by 3 list of triangle indices into V
// Outputs:
//   VAO  identifier of compiled vertex array object.
inline void mesh_to_vao(
  const GLuint & id,
  const Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> & V,
  const Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> & F,
  const Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> & N,
  const Eigen::Matrix< float,Eigen::Dynamic,2,Eigen::RowMajor> & TC,
  GLuint & VAO);

// Implementation

inline void mesh_to_vao(
  const GLuint & id,
  const Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> & V,
  const Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> & F,
  const Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> & N,
  const Eigen::Matrix< float,Eigen::Dynamic,2,Eigen::RowMajor> & TC,
  GLuint & VAO)
{
    glGenVertexArrays(1, &VAO);
    // Generate and attach buffers to vertex array
    GLuint VBO, EBO, NBO, TCBO;//, FBO;
    //GLuint depth_render_buf, color_render_buf;
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glGenBuffers(1, &NBO);
    glGenBuffers(1, &TCBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*V.size(), V.data(), GL_STATIC_DRAW);
    
/*
    glGenVertexArrays(1, &planeVAO);
    glGenBuffers(1, &planeVBO);
    glBindVertexArray(planeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, planeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(planeVertices), planeVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glBindVertexArray(0);
*/

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*F.size(), F.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

    // glBindBuffer(GL_ARRAY_BUFFER, NBO);
    // glBufferData(GL_ARRAY_BUFFER, sizeof(GLuint)*N.size(), N.data(), GL_STATIC_DRAW);
    // glEnableVertexAttribArray(1);
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

    bool not_empty = (TC.array() != 0.0).any();
    if(not_empty)
    {
      //std::cout << "not empty!" << std::endl;
      glBindBuffer(GL_ARRAY_BUFFER, TCBO);
      glBufferData(GL_ARRAY_BUFFER, sizeof(GLuint)*V.size(), TC.data(), GL_STATIC_DRAW);
      glEnableVertexAttribArray(1);
      glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), (GLvoid*)0);
    }

    glBindBuffer(GL_ARRAY_BUFFER, 0); 
    glBindVertexArray(0);

    // Generate and attach frame buffers
    // glGenFramebuffers(1, &FBO);

    // glGenRenderbuffers(1, &color_render_buf);
    // glBindRenderbuffer(GL_RENDERBUFFER, color_render_buf);
    // glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, w, h);
    // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, color_render_buf);

    // glGenRenderbuffers(1,&depth_render_buf);
    // glBindRenderbuffer(GL_RENDERBUFFER, depth_render_buf);
    // glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, w, h);
    // glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depth_render_buf);

    // glBindFramebuffer(0, FBO);

}

#endif