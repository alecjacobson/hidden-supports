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
  const int & height,
  const std::string type);

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
  const int & height,
  const std::string type)
{
    glGenFramebuffers(1, &fbo);

    // Create the depth buffer
    glGenTextures(1, &shadow_map);
    glActiveTexture(id);
    glBindTexture(GL_TEXTURE_2D, shadow_map);
    if(type == "depth")
      glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
    else if(type == "color" or type == "none")
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width, height, 0, GL_RED, GL_FLOAT, NULL);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);

    igl::opengl::report_gl_error("gen textures\n");


    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    igl::opengl::report_gl_error("bind fbo\n");

    if(type == "depth")
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, shadow_map, 0);
    else if(type == "color")
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, shadow_map, 0);
    else if(type == "none")
      std::cout << "none" << std::endl;
    igl::opengl::report_gl_error("attachments\n");

    // glDrawBuffer(GL_NONE);
    // glReadBuffer(GL_NONE);

    if(type == "depth" || type == "color")
    {
      bool status = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
      if(!status)
          std::cout << "Could not initialize FBO :(" << std::endl;
      // else
          // std::cout << "FBO ready!" << std::endl;

    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);


}

inline void bind_map_for_writing(
  GLuint & fbo)
{
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    // glReadBuffer(GL_NONE);

}

inline void bind_map_for_reading(
  GLuint & shadow_map,
  GLenum TextureUnit)
{
    glActiveTexture(TextureUnit);
    glBindTexture(GL_TEXTURE_2D, shadow_map);
}

#endif