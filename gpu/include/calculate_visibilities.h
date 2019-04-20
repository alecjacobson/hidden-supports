#ifndef CALCULATE_VISIBILITIES_H
#define CALCULATE_VISIBILITIES_H

// // make sure the modern opengl headers are included before any others
// #include "gl.h"

// // #include <OpenGL/gl3.h>
// // #define __gl_h_
// #define GLFW_INCLUDE_GLU
// #include <GLFW/glfw3.h>
#include "helpers.h"

#include "perspective_projection.h"
#include "orthographic_projection.h"
#include "print_opengl_info.h"
#include "report_gl_error.h"
#include "create_shader_program_from_files.h"
#include "last_modification_time.h"
#include "mesh_to_vao.h"
#include "shadow_buffer.h"
#include "make_light_matrix.h"
#include "read_json.h"

inline int calculate_visibilities(
  const Eigen::Vector3i &side,  // size of scene (voxel grid)
  const Eigen::MatrixXf &views, // viewpoint origins
  const std::vector<std::string> paths,
  const scene &sc,
  Eigen::VectorXi &visibilities
);

inline int calculate_visibilities(
  const Eigen::Vector3i &side,  // size of scene (voxel grid)
  const Eigen::MatrixXf &views, // viewpoint origins
  const std::vector<std::string> paths,
  const scene &sc,
  Eigen::VectorXi &visibilities
)
{
  GLuint prog_id=0;
  GLuint q_prog_id=0;
  GLuint render_prog_id=0;
  GLuint VAO,Q_VAO;
  GLuint FBO,FBO_render_even,FBO_render_odd,FBO_large;
  GLuint shadow_map,visibility_map_even,visibility_map_odd;

  Eigen::Affine3f light_view = 
    Eigen::Affine3f::Identity() * 
    Eigen::Translation3f(Eigen::Vector3f(0,0,-1));

  std::vector<std::string> vertex_shader_paths;
  std::vector<std::string> fragment_shader_paths;
  std::vector<std::string> q_vertex_shader_paths;
  std::vector<std::string> q_fragment_shader_paths;
  std::vector<std::string> render_vertex_shader_paths;
  std::vector<std::string> render_fragment_shader_paths;

  int w = side(0);
  int h = side(1);
  int d = side(2);

  float near = 0.1;
  float far = 1000;
  float z_slice = 0.0;
  float light_top = tan((40./2.)*M_PI/180.)*near;
  float top = 0.5;
  float light_right = light_top * (double)w/(double)h;
  float right = top * (double)w/(double)h;

  float max_z = 0.0;
  float min_z = 0.0;
  float z_range = 0.0;
  float step = 0.0;

  Eigen::Affine3f view = 
    Eigen::Affine3f::Identity() * 
    Eigen::Translation3f(Eigen::Vector3f(0,0,-1));

  Eigen::Matrix4f proj = Eigen::Matrix4f::Identity();
  Eigen::Matrix4f light_proj = Eigen::Matrix4f::Identity();
  Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
  Eigen::Matrix4f q_model = Eigen::Matrix4f::Identity();

  if(!glfwInit())
  {
     std::cerr<<"Could not initialize glfw"<<std::endl;
     return 0;
  }
  const auto & error = [] (int error, const char* description)
  {
    std::cerr<<description<<std::endl;
  };

  glfwSetErrorCallback(error);
  glfwWindowHint(GLFW_SAMPLES, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

  std::cout << "width and height " << w << " " << h << std::endl;
  GLFWwindow* window = glfwCreateWindow(w, h, "visibility", NULL, NULL);

  if(!window)
  {
    glfwTerminate();
    std::cerr<<"Could not create glfw window"<<std::endl;
    return 0;
  }

  glfwSetWindowPos(window,0,0);
  glfwMakeContextCurrent(window);

  // Load OpenGL and its extensions
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
  {
    std::cerr<<"Failed to load OpenGL and its extensions"<<std::endl;
    return 0;
  }
  print_opengl_info(window);
  igl::opengl::report_gl_error("init");

  glfwSetInputMode(window,GLFW_CURSOR,GLFW_CURSOR_NORMAL);

  GLint v_dims[2];
  GLint t_dims[1];
  glGetIntegerv(GL_MAX_VIEWPORT_DIMS, &v_dims[0]);
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &t_dims[0]);
  std::cout << v_dims[0] << " " << v_dims[1] << std::endl;
  std::cout << t_dims[0] << std::endl;

  Eigen::Vector3f centroid;
  igl::centroid(sc.V,sc.F,centroid);

  int t_w = w*2;
  int t_h = h*2;
  // ratio = (float)std::max(w,h) / (float)std::min(w,h);
  float ratio = (float)w / (float)h;
  std::cout << side << std::endl;

  max_z = sc.V.col(2).maxCoeff();
  min_z = sc.V.col(2).minCoeff();
  z_range = std::abs(max_z - min_z);
  int number_of_slices = side(2);
  step = z_range / number_of_slices;

  // Projection and modelview matrices
  perspective(-light_right, light_right, -light_top, light_top, near, far, light_proj);
  // Orthographic projection matrix
  orthographic(-right, right, -top, top, near, far, proj);

  Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> Q_V;
  Eigen::Matrix< float,Eigen::Dynamic,2,Eigen::RowMajor> Q_TC;
  Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> Q_N;

  Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> Q_F;
  Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> Q_FTC;
  Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> Q_FN;

  igl::readOBJ("../data/u-quad.obj", Q_V, Q_TC, Q_N, Q_F, Q_FTC, Q_FN);
  Q_V.rowwise() -= Q_V.colwise().mean();
  Q_V /= (Q_V.colwise().maxCoeff()-Q_V.colwise().minCoeff()).maxCoeff();

  // Close the window if user presses ESC or CTRL+C
  glfwSetKeyCallback(
    window,
    [](GLFWwindow* window, int key, int scancode, int action, int mods)
    {
      if(key == 256 || (key == 67 && (mods & GLFW_MOD_CONTROL)))
      {
        glfwSetWindowShouldClose(window,true);
      }
  });


  glEnable(GL_DEPTH_TEST);

  // Force compilation on first iteration through loop
  double time_of_last_shader_compilation = 0;
  double time_of_last_json_load = 0;
  const auto any_changed = 
    [](
        const std::vector<std::string> &paths,
        const double time_of_last_shader_compilation
        )->bool
  {
    for(const auto & path : paths)
    {
      if(last_modification_time(path) > time_of_last_shader_compilation)
      {
        std::cout<<path<<" has changed since last compilation attempt."<<std::endl;
        return true;
      }
    }
    return false;
  };

  float start_time = igl::get_seconds();

  z_slice = min_z;
  // z_slice = min_z + (number_of_slices/2)*step;
  // max_z = z_slice + step;

  int count = 0;
  // int copy_count = 0;
  // int which_texture = 0;
  // int tex_dims = t_dims[0];
  // GLint y_offset = 0;

  Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_values;
  visibility_values.resize(w*h*number_of_slices,1);

  // Main display routine
  const auto & tictoc = []()
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };

  while(z_slice < max_z-step)
  {
    std::cout << "Z SLICE NUMBER: " << count << std::endl;
    std::cout << "Z SLICE DEPTH: " << z_slice << std::endl;
    std::cout << "-----" << std::endl;



    if(any_changed(paths,time_of_last_json_load))
    {
      std::cout<<"-----------------------------------------------"<<std::endl;
      time_of_last_json_load = igl::get_seconds();
      if(!read_json(paths[0],
            vertex_shader_paths,
            fragment_shader_paths))
      {
        std::cerr<<"Failed to read "<<paths[0]<<std::endl;
      }
      if(!read_json(paths[1],
            q_vertex_shader_paths,
            q_fragment_shader_paths))
      {
        std::cerr<<"Failed to read "<<paths[1]<<std::endl;
      }
      if(!read_json(paths[2],
            render_vertex_shader_paths,
            render_fragment_shader_paths))
      {
        std::cerr<<"Failed to read "<<paths[2]<<std::endl;
      }
      // force reload of shaders
      time_of_last_shader_compilation = 0;
    }
    if(
      any_changed(vertex_shader_paths         ,time_of_last_shader_compilation) ||
      any_changed(fragment_shader_paths       ,time_of_last_shader_compilation) ||
      any_changed(q_vertex_shader_paths         ,time_of_last_shader_compilation)||
      any_changed(q_fragment_shader_paths       ,time_of_last_shader_compilation)||
      any_changed(render_vertex_shader_paths    ,time_of_last_shader_compilation)||
      any_changed(render_fragment_shader_paths  ,time_of_last_shader_compilation)) 
    {
      std::cout<<"-----------------------------------------------"<<std::endl;
      // remember the time we tried to compile
      time_of_last_shader_compilation = igl::get_seconds();
      if(
          !create_shader_program_from_files(
            vertex_shader_paths,
            fragment_shader_paths,
            prog_id))
      {
        // Force null shader to visually indicate failure
        glDeleteProgram(prog_id);
        prog_id = 0;
        std::cout<<"-----------------------------------------------"<<std::endl;
      }
      if(
          !create_shader_program_from_files(
            q_vertex_shader_paths,
            q_fragment_shader_paths,
            q_prog_id))
      {
        // Force null shader to visually indicate failure
        glDeleteProgram(q_prog_id);
        q_prog_id = 0;
        std::cout<<"-----------------------------------------------"<<std::endl;
      }
      if(
          !create_shader_program_from_files(
            render_vertex_shader_paths,
            render_fragment_shader_paths,
            render_prog_id))
      {
        // Force null shader to visually indicate failure
        glDeleteProgram(render_prog_id);
        render_prog_id = 0;
        std::cout<<"-----------------------------------------------"<<std::endl;
      }
    }
    igl::opengl::report_gl_error("loaded shaders\n");

    init_shadow_buffer(shadow_map, FBO, GL_TEXTURE0, t_w, t_h, "depth");
    igl::opengl::report_gl_error("init shadow buffer\n");

    init_shadow_buffer(visibility_map_odd, FBO_render_odd, GL_TEXTURE1, w, h, "color");
    igl::opengl::report_gl_error("init shadow buffer\n");

    init_shadow_buffer(visibility_map_even, FBO_render_even, GL_TEXTURE1, w, h, "color");
    igl::opengl::report_gl_error("init shadow buffer\n");

    for(int v = 0; v < views.rows(); v++)
    {
      Eigen::Vector3f viewpoint = views.row(v);

      make_light_matrix(viewpoint,centroid,light_view);

      // select program and attach uniforms
      glUseProgram(prog_id);

      GLint light_proj_loc = glGetUniformLocation(prog_id,"light_proj");
      glUniformMatrix4fv(light_proj_loc,1,GL_FALSE,light_proj.data());
      GLint proj_loc = glGetUniformLocation(prog_id,"proj");
      glUniformMatrix4fv(proj_loc,1,GL_FALSE,proj.data());
      GLint model_loc = glGetUniformLocation(prog_id,"model");
      glUniformMatrix4fv(model_loc,1,GL_FALSE,model.data());
      GLint light_view_loc = glGetUniformLocation(prog_id,"light_view");
      // glUniformMatrix4fv(light_view_loc,num_views,GL_FALSE,light_views[0]);
      glUniformMatrix4fv(light_view_loc,1,GL_FALSE,light_view.data());
      GLint p_view_loc = glGetUniformLocation(prog_id,"view");
      glUniformMatrix4fv(p_view_loc,1,GL_FALSE,view.data());

      igl::opengl::report_gl_error("uniforms using prog_id\n");

      //////////

      bind_map_for_writing(FBO);

      glClear(GL_DEPTH_BUFFER_BIT);

      mesh_to_vao(prog_id, sc.V, sc.F, sc.N, sc.TC, VAO);
      igl::opengl::report_gl_error("bind vao 1\n");
      
      glBindVertexArray(VAO);

      glViewport(0, 0, t_w, t_h);

      // draw elements to texture
      glDrawElements(GL_TRIANGLES, sc.F.size(), GL_UNSIGNED_INT, 0);
      igl::opengl::report_gl_error("draw elements to texture\n");

      // Eigen::Matrix< GLubyte,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> shadow;
      // shadow.resize(t_w*t_h,1);

      // glReadPixels(0, 0, t_w, t_h, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, shadow.data());
      // igl::writeDMAT("texture"+std::to_string(v)+".dmat", shadow,true);

      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      igl::opengl::report_gl_error("unbind fbo\n");

      glViewport(0, 0, w, h);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      ///////////
      q_model << ratio, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, z_slice,
                  0, 0, 0, 1;

      // select program and attach uniforms
      glUseProgram(q_prog_id);
      GLint q_texture_loc = glGetUniformLocation(q_prog_id,"shadow_map");
      glUniform1i(q_texture_loc, 0);
      GLint q_near_loc = glGetUniformLocation(q_prog_id,"near_plane");
      glUniform1f(q_near_loc, near);
      GLint q_far_loc = glGetUniformLocation(q_prog_id,"far_plane");
      glUniform1f(q_far_loc, far);
      GLint q_proj_loc = glGetUniformLocation(q_prog_id,"proj");
      glUniformMatrix4fv(q_proj_loc,1,GL_FALSE,proj.data());
      GLint q_light_proj_loc = glGetUniformLocation(q_prog_id,"light_proj");
      glUniformMatrix4fv(q_light_proj_loc,1,GL_FALSE,light_proj.data());
      GLint q_model_loc = glGetUniformLocation(q_prog_id,"model");
      glUniformMatrix4fv(q_model_loc,1,GL_FALSE,model.data());
      GLint q_q_model_loc = glGetUniformLocation(q_prog_id,"q_model");
      glUniformMatrix4fv(q_q_model_loc,1,GL_FALSE,q_model.data());
      GLint q_view_loc = glGetUniformLocation(q_prog_id,"view");
      glUniformMatrix4fv(q_view_loc,1,GL_FALSE,view.data());
      GLint q_light_view_loc = glGetUniformLocation(q_prog_id,"light_view");
      glUniformMatrix4fv(q_light_view_loc,1,GL_FALSE,light_view.data());

      igl::opengl::report_gl_error("uniforms using q_prog_id\n");
      
      if(v == 0)
      {
        // first one
        bind_map_for_writing(FBO_render_even);
        GLint q_index_loc = glGetUniformLocation(q_prog_id,"index");
        glUniform1i(q_index_loc, v);
      }
      else if(v % 2 == 0)
      {
        GLint q_vtexture_loc = glGetUniformLocation(q_prog_id,"visibility_map");
        glUniform1i(q_vtexture_loc, 1);
        bind_map_for_writing(FBO_render_even);
        bind_map_for_reading(visibility_map_odd, GL_TEXTURE1);
        GLint q_index_loc = glGetUniformLocation(q_prog_id,"index");
        glUniform1i(q_index_loc, v);
      }
      else if(v % 2 != 0)
      {
        GLint q_vtexture_loc = glGetUniformLocation(q_prog_id,"visibility_map");
        glUniform1i(q_vtexture_loc, 1);
        bind_map_for_writing(FBO_render_odd);
        bind_map_for_reading(visibility_map_even, GL_TEXTURE1);
        GLint q_index_loc = glGetUniformLocation(q_prog_id,"index");
        glUniform1i(q_index_loc, v);
      }
      igl::opengl::report_gl_error("binding maps for reading/writing\n");

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
      mesh_to_vao(q_prog_id, Q_V, Q_F, Q_N, Q_TC,Q_VAO);

      bind_map_for_reading(shadow_map, GL_TEXTURE0);

      glBindVertexArray(Q_VAO);
      igl::opengl::report_gl_error("bind vao quad\n");

      glDrawElements(GL_TRIANGLES, Q_F.size(), GL_UNSIGNED_INT, 0);
      igl::opengl::report_gl_error("draw elements to offscreen buffer\n");
      
      // Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slice;
      // visibility_slice.resize(w*h,1);
      // glReadPixels(0, 0, w, h, GL_RED, GL_FLOAT, visibility_slice.data());
      // igl::writeDMAT("slice"+std::to_string(v)+".dmat", visibility_slice, true);

      glBindVertexArray(0);

      ///////////
      glfwSwapBuffers(window);
      
    }

    igl::opengl::report_gl_error("delete\n");

    // end of view for loop

    Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slice;
    visibility_slice.resize(w*h,1);

    glReadPixels(0, 0, w, h, GL_RED, GL_FLOAT, visibility_slice.data());
    // igl::writeDMAT("slice"+std::to_string(count)+".dmat", visibility_slice, true);

    if(count <= side(2))
    {
      visibility_values.block(count*w*h, 0, w*h, 1) = visibility_slice;
    }

    z_slice += step;
    count++;

    glDeleteTextures(1,&shadow_map);
    glDeleteTextures(1,&visibility_map_even);
    glDeleteTextures(1,&visibility_map_odd);

    glDeleteFramebuffers(1,&FBO);
    glDeleteFramebuffers(1,&FBO_render_odd);
    glDeleteFramebuffers(1,&FBO_render_even);
    
  }


  std::cout << "size of matrix: " << visibility_values.rows() << std::endl;
  Eigen::MatrixXi vv = visibility_values.cast <int> ();
  Eigen::VectorXi S(Eigen::Map<Eigen::VectorXi>(vv.data(), vv.cols()*vv.rows()));  
  std::cout << "size of vector: " << S.rows() << std::endl;
  visibilities = S;
  // igl::writeDMAT("../results/visibilities.dmat", S, true);

  glDeleteVertexArrays(1,&VAO);
  glDeleteFramebuffers(1,&FBO_large);


  glfwDestroyWindow(window);
  glfwTerminate();

  return 1;
}

#endif