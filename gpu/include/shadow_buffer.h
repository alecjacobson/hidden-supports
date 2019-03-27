#ifndef SHADOW_BUFFER_H
#define SHADOW_BUFFER_H
#include "gl.h"
#include <Eigen/Core>
#include <iostream>

inline void init_shadow_buffer(
  GLuint & shadow_map,
  GLuint & fbo,
  const GLuint & id,
  const int & width,
  const int & height);

inline void bind_map_for_reading(
  GLuint & shadow_map,
  GLenum TextureUnit);

inline void bind_map_for_writing(
  GLuint & fbo);

// Implementation

inline void init_shadow_buffer(
  GLuint & shadow_map,
  GLuint & fbo,
  const GLuint & id,
  const int & width,
  const int & height)
{
    glGenFramebuffers(1, &fbo);

    // Create the depth buffer
    glGenTextures(1, &shadow_map);
    glBindTexture(GL_TEXTURE_2D, shadow_map);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);

    igl::opengl::report_gl_error("gen textures\n");


    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    igl::opengl::report_gl_error("bind fbo\n");

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadow_map, 0);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);

    bool status = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    if(!status)
        std::cout << "Could not initialise FBO" << std::endl;
    else
        std::cout << "FBO ready!" << std::endl;

}

inline void bind_map_for_writing(
  GLuint & fbo)
{
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    //glReadBuffer(GL_NONE);

}

inline void bind_map_for_reading(
  GLuint & shadow_map,
  GLenum TextureUnit)
{
    glActiveTexture(TextureUnit);
    glBindTexture(GL_TEXTURE_2D, shadow_map);
}

#endif