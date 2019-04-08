#ifndef RENDER_BUFFER_H
#define RENDER_BUFFER_H
#include "gl.h"
#include <Eigen/Core>
#include <iostream>

inline void init_render_buffer(
  GLuint & color_map,
  GLuint & fbo,
  const GLuint & id,
  const int & width,
  const int & height,
  const std::string type);

inline void init_render_buffer(
  GLuint & color_map,
  GLuint & fbo,
  const GLuint & id,
  const int & width,
  const int & height,
  const std::string type)
  {
  render buffer as color buffer
  glGenFramebuffers(1, &fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo);
  glGenFramebuffers(1, &fbo);
  glGenRenderbuffers(1, &color_map);

  glBindRenderbuffer(GL_RENDERBUFFER, color_map);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_RED, w, h);
  glBindFramebuffer(GL_FRAMEBUFFER, fbo_render);

  attach render buffer to the fbo as color buffer
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, color_map);
  bool status = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
  if(!status)
      std::cout << "Could not initialise FBO" << std::endl;
  else
      std::cout << "color FBO ready!" << std::endl;
  }