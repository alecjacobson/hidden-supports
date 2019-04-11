// make sure the modern opengl headers are included before any others
#include "gl.h"

// #include <OpenGL/gl3.h>
// #define __gl_h_
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include <igl/frustum.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <igl/get_seconds.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/remove_duplicate_vertices.h>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/centroid.h>

#include <string>
#include <chrono>
#include <thread>
#include <iostream>
#include <climits>

#include "read_json.h"
#include "print_opengl_info.h"
#include "report_gl_error.h"
#include "create_shader_program_from_files.h"
#include "last_modification_time.h"
#include "mesh_to_vao.h"
#include "shadow_buffer.h"
#include "perspective_projection.h"
#include "orthographic_projection.h"
#include "look_at.h"
#include "voxelize.h"
#include "generate_views.h"

// int w=512,h=301;
int w=200,h=139;
int d=200;
int t_w,t_h;
float ratio = 0.0;

// max viewport size 16384
// 16384^(1/1.5) ~= 645
// 24^2 = 576
//double highdpi=1;
GLuint prog_id=0;
GLuint q_prog_id=0;
GLuint render_prog_id=0;
GLuint VAO,Q_VAO;
GLuint FBO,FBO_render_even,FBO_render_odd,FBO_large;
GLuint shadow_map,visibility_map_even,visibility_map_odd;

bool wire_frame = false;
bool mouse_down = false;

// Eigen::Matrix4f view = Eigen::Matrix4f::Identity();
Eigen::Affine3f view = 
  Eigen::Affine3f::Identity() * 
  Eigen::Translation3f(Eigen::Vector3f(0,0,-1));

Eigen::Matrix4f proj = Eigen::Matrix4f::Identity();
Eigen::Matrix4f light_proj = Eigen::Matrix4f::Identity();
Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
Eigen::Matrix4f q_model = Eigen::Matrix4f::Identity();

Eigen::Affine3f light_view = 
  Eigen::Affine3f::Identity() * 
  Eigen::Translation3f(Eigen::Vector3f(0,0,-1));

// Eigen::Matrix4f light_view = Eigen::Matrix4f::Identity();

// Eigen::Vector4f eye = Eigen::Vector4f(4,4,5,1);
// Eigen::Vector4f at = Eigen::Vector4f(4,4,-4.5,1);
// Eigen::Vector4f up = Eigen::Vector4f(0,1,0,0);

float near = 0.1;
// float near = 1.0;
float far = 1000;
float z_slice = 0.0;
float light_top = tan((60./2.)*M_PI/180.)*near;
float top = 0.5;
float light_right = light_top * (double)::w/(double)::h;
float right = top * (double)::w/(double)::h;

float max_z = 0.0;
float min_z = 0.0;
float z_range = 0.0;
float step = 0.0;

// Mesh data: RowMajor is important to directly use in OpenGL
Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> V;
Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> F;
Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> N;
Eigen::Matrix< float,Eigen::Dynamic,2,Eigen::RowMajor> TC;

Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> Q_V;
Eigen::Matrix< float,Eigen::Dynamic,2,Eigen::RowMajor> Q_TC;
Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> Q_N;

Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> Q_F;
Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> Q_FTC;
Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> Q_FN;


int main(int argc, char * argv[])
{
  Eigen::setNbThreads(1);

  std::vector<std::string> vertex_shader_paths;
  std::vector<std::string> fragment_shader_paths;
  std::vector<std::string> q_vertex_shader_paths;
  std::vector<std::string> q_fragment_shader_paths;
  std::vector<std::string> render_vertex_shader_paths;
  std::vector<std::string> render_fragment_shader_paths;

  if(!glfwInit())
  {
     std::cerr<<"Could not initialize glfw"<<std::endl;
     return EXIT_FAILURE;
  }
  const auto & error = [] (int error, const char* description)
  {
    std::cerr<<description<<std::endl;
  };

  glfwSetErrorCallback(error);
  glfwWindowHint(GLFW_SAMPLES, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  //glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);

  GLFWwindow* window = glfwCreateWindow(w, h, "visibility", NULL, NULL);

  if(!window)
  {
    glfwTerminate();
    std::cerr<<"Could not create glfw window"<<std::endl;
    return EXIT_FAILURE;
  }

  glfwSetWindowPos(window,0,0);
  glfwMakeContextCurrent(window);

  // Load OpenGL and its extensions
  if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
  {
    std::cerr<<"Failed to load OpenGL and its extensions"<<std::endl;
    return EXIT_FAILURE;
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

  // Read input scene from file
  igl::readSTL(argv[4], V, F, N);

  Eigen::RowVector3f translation = V.colwise().mean();
  V.rowwise() -= translation;

  float scale_factor = (V.colwise().maxCoeff()-V.colwise().minCoeff()).maxCoeff();
  V /= scale_factor;

  Eigen::Vector3f centroid;
  igl::centroid(V,F,centroid);
  // look_at(Eigen::Vector3f(0,0,-1), Eigen::Vector3f(0,1,0), centroid, view);

  // make voxel grid
  Eigen::MatrixXf GV;
  Eigen::Vector3i side;
  Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor>  F_int = F.cast <int> ();
  std::cout << "casted to int" << std::endl;
  surround_scene_in_grid(d, V, F_int, side, GV);
  w = side(0);
  h = side(1);
  d = side(2);

  t_w = w*2;
  t_h = h*2;
  ratio = (float)std::max(w,h) / (float)std::min(w,h);
  // ratio = (float)w / (float)h;
  std::cout << "made voxel grid" << std::endl;
  std::cout << side << std::endl;
  // std::cout << ratio << std::endl;

  max_z = V.col(2).maxCoeff();
  min_z = V.col(2).minCoeff();
  // std::cout << max_z << std::endl << min_z << std::endl;
  z_range = std::abs(max_z - min_z);
  int number_of_slices = side(2)+1;
  step = z_range / number_of_slices;

  // Projection and modelview matrices
  perspective(-light_right, light_right, -light_top, light_top, near, far, light_proj);
  // Orthographic projection matrix
  orthographic(-right, right, -top, top, near, far, proj);
  // perspective(-light_right, light_right, -light_top, light_top, near, far, proj);

  // get view rays
  float num_views = 2.0;
  Eigen::Vector3f bottom_left = V.colwise().minCoeff();
	Eigen::Vector3f top_right = V.colwise().maxCoeff();
	Eigen::MatrixXf views;
	// generate_views(bottom_left, top_right, num_views, views);
  views.resize(2,3);
  views.row(0) = Eigen::Vector3f(0,0,-1);
  views.row(1) = Eigen::Vector3f(-0.1,0,-1);
	std::cout << "generated viewing rays" << std::endl;

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


  float fl_num_textures = (float)(w*h) / (float)(t_dims[0]);
  std::cout << "FLOAT NUM TEXTURES " << fl_num_textures << std::endl;

  int num_textures = std::ceil(fl_num_textures);
  std::cout << "NUM TEXTURES " << num_textures << std::endl;

  float fl_max_slices_per_texture = (number_of_slices) / num_textures;
  int max_slices_per_texture = std::floor(fl_max_slices_per_texture);

  GLuint large_visibilities[num_textures];

  if(num_textures > 1)
  {
    for(int i = 0; i < num_textures; i++)
    {
      init_shadow_buffer(large_visibilities[i], FBO_large, render_prog_id, 
          w, h*max_slices_per_texture, "none");
      std::cout << "SIZE " << w << " by " << h*max_slices_per_texture << std::endl;
    }
  }
  else
  {
    init_shadow_buffer(large_visibilities[0], FBO_large, render_prog_id, 
        w, h*(number_of_slices), "none");
    igl::opengl::report_gl_error("init large shadow buffer\n");
    std::cout << "SIZE " << w << " by " << h*(number_of_slices) << std::endl;
  }

  igl::opengl::report_gl_error("made large textures\n");

  float start_time = igl::get_seconds();

  z_slice = min_z;
  // z_slice = min_z + (number_of_slices/2)*step;
  // max_z = z_slice + step;

  int count = 0;
  int copy_count = 0;
  int which_texture = 0;
  int tex_dims = t_dims[0];
  GLint y_offset = 0;

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

  while(z_slice < max_z)
  {
    std::cout << "Z SLICE NUMBER: " << count << std::endl;
    std::cout << "Z SLICE DEPTH: " << z_slice << std::endl;
    std::cout << "-----" << std::endl;



    if(any_changed({argv[1]},time_of_last_json_load))
    {
      std::cout<<"-----------------------------------------------"<<std::endl;
      time_of_last_json_load = igl::get_seconds();
      if(!read_json(argv[1],
            vertex_shader_paths,
            fragment_shader_paths))
      {
        std::cerr<<"Failed to read "<<argv[1]<<std::endl;
      }
      if(!read_json(argv[2],
            q_vertex_shader_paths,
            q_fragment_shader_paths))
      {
        std::cerr<<"Failed to read "<<argv[2]<<std::endl;
      }
      if(!read_json(argv[3],
            render_vertex_shader_paths,
            render_fragment_shader_paths))
      {
        std::cerr<<"Failed to read "<<argv[3]<<std::endl;
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

      init_shadow_buffer(shadow_map, FBO, q_prog_id, t_w, t_h, "depth");
      igl::opengl::report_gl_error("init shadow buffer\n");

      init_shadow_buffer(visibility_map_odd, FBO_render_odd, q_prog_id, w, h, "color");
      igl::opengl::report_gl_error("init shadow buffer\n");

      init_shadow_buffer(visibility_map_even, FBO_render_even, q_prog_id, w, h, "color");
      igl::opengl::report_gl_error("init shadow buffer\n");

      for(int v = 0; v < views.rows(); v++)
      {
        

        std::cout << "view tic " << v << " " << tictoc() << std::endl;

        Eigen::Vector3f viewpoint = views.row(v);
        Eigen::Vector3f l = centroid - viewpoint;
        l.normalize();
        Eigen::Vector3f n = -(viewpoint - Eigen::Vector3f(0,0,1));
        n.normalize();
        float cos_spin_angle = n.dot(l);
        // if x value is negative, then rotate positive angle
        float spin_angle = (viewpoint(0) < 0) ? acos(cos_spin_angle) : -1*acos(cos_spin_angle);

        light_view = Eigen::Affine3f::Identity() * 
          Eigen::Translation3f(views.row(v));
        light_view.rotate(
          Eigen::AngleAxisf(
            spin_angle,
            Eigen::Vector3f(0,1,0)));
        // std::cout << "view " << views.row(v) << std::endl; 
        // look_at(views.row(v), Eigen::Vector3f(0,1,0), centroid, light_view);

        // clear screen and set viewport
        // glClearColor(1.0,1.0,1.0,0.);
        // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // igl::opengl::report_gl_error("cleared screen\n");


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

        mesh_to_vao(prog_id, V, F, N, TC,VAO);
        igl::opengl::report_gl_error("bind vao 1\n");
        
        glBindVertexArray(VAO);

        glViewport(0, 0, t_w, t_h);

        // draw elements to texture
        glDrawElements(GL_TRIANGLES, F.size(), GL_UNSIGNED_INT, 0);
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

        // glFlush();
        // glfwSwapBuffers(window);

        glBindVertexArray(0);

        ///////////

        // glFlush();
        glfwSwapBuffers(window);
        // igl::opengl::report_gl_error("flush\n");

      
        
      }

      igl::opengl::report_gl_error("delete\n");

      // end of view for loop

      Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slice;
      visibility_slice.resize(w*h,1);


      // if (copy_count >= max_slices_per_texture)
      // { 
      //   which_texture++; 
      //   tex_dims = t_dims[0] * (which_texture+1);
      //   y_offset = 0;
      //   copy_count = 0;
      // }
      
      // std::cout << "WHICH TEXTURE? " << which_texture << std::endl;
      // // std::cout << y_offset << std::endl;

      // bind_map_for_reading(large_visibilities[which_texture], GL_TEXTURE2);
      // glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, y_offset, 0, 0, w, h);
      // // std::cout << "copy texture tic " << tic << std::endl;
      // igl::opengl::report_gl_error("copy texture\n");

      // y_offset += h;


      glReadPixels(0, 0, w, h, GL_RED, GL_FLOAT, visibility_slice.data());
      igl::writeDMAT("slice"+std::to_string(count)+".dmat", visibility_slice, true);

      if(count <= side(2))
      {
        visibility_values.block(count*w*h, 0, w*h, 1) = visibility_slice;
      }

      z_slice += step;
      count++;
      copy_count++;

      glDeleteTextures(1,&shadow_map);
      glDeleteTextures(1,&visibility_map_even);
      glDeleteTextures(1,&visibility_map_odd);

      glDeleteFramebuffers(1,&FBO);
      glDeleteFramebuffers(1,&FBO_render_odd);
      glDeleteFramebuffers(1,&FBO_render_even);
      
    }

    // Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slices;
    // visibility_slices.resize((number_of_slices)*w*h,1);
    // for(int i = 0; i < num_textures; i++)
    // {
    //     bind_map_for_reading(large_visibilities[i],GL_TEXTURE2);

    //     Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slices;
    //     visibility_slices.resize(max_slices_per_texture*w*h,1);

    //     glGetTexImage(GL_TEXTURE_2D, 0, GL_RED, GL_FLOAT, visibility_slices.data());
    //     igl::opengl::report_gl_error("get tex\n");
    //     // glClearTexImage(large_visibilities[i], 0, GL_RED, GL_FLOAT, NULL);
    //     // igl::opengl::report_gl_error("clear tex\n");

    //     // std::cout << "BLOCK " << "0, " << i*max_slices_per_texture*w*h << ", " 
    //     //           << max_slices_per_texture*w*h << ", 1" << 
    //     //           " to 0, 0, " <<max_slices_per_texture*w*h << ", 1" << std::endl;

    //     visibility_values.block(
    //               i*max_slices_per_texture*w*h, 0,
    //               max_slices_per_texture*w*h, 1)
    //               = visibility_slices.block(0,0,max_slices_per_texture*w*h, 1);
        
    //     // igl::writeDMAT("slices"+std::to_string(i)+".dmat", visibility_slices, true);
    // }

    glDeleteVertexArrays(1,&VAO);
    glDeleteFramebuffers(1,&FBO_large);

    std::cout << "size of vector: " << visibility_values.rows() << std::endl;
    Eigen::VectorXf S(Eigen::Map<Eigen::VectorXf>(visibility_values.data(), 
            visibility_values.cols()*visibility_values.rows()));      
    // all invisible/geometry voxels have 1 or floats, visible voxels have 0
    S /= num_views;
    igl::writeDMAT("visibilities.dmat", S, true);
    // std::cout << "filled voxels from main " << (S.array() == 1).count() << std::endl;

    Eigen::MatrixXf V_voxels;
    Eigen::MatrixXi F_voxels;
    Eigen::MatrixXf V_voxels_mc;
    Eigen::MatrixXi F_voxels_mc;
    std::cout << side << std::endl;

    double isovalue = 1.0;
    // voxelize
    make_voxels_from_visibility(S, GV, side, isovalue, V_voxels, F_voxels);
    std::cout << V_voxels.rows() << ", " << V_voxels.cols() << std::endl;
    std::cout << F_voxels.rows() << ", " << F_voxels.cols() << std::endl;
    // write obj
    /////////
    Eigen::Vector3f m = V.colwise().maxCoeff();
    V *= scale_factor;
    V.rowwise() += translation;

    V_voxels *= scale_factor;
    V_voxels.rowwise() += translation;

    GV *= scale_factor;
    GV.rowwise() += translation;
    /////////////////
    Eigen::MatrixXf NV;
    Eigen::MatrixXi NF;
    Eigen::MatrixXi IM;
    Eigen::MatrixXi MI;
    igl::remove_duplicate_vertices(V_voxels, F_voxels, 1e-8, NV, IM, MI, NF);

    igl::writeOBJ("output_voxels.obj", NV, NF);
    // igl::writeOBJ("output_voxels.obj", V_voxels, F_voxels);

    // Eigen::VectorXf sum = V_voxels.rowwise().sum();
    // std::cout << (sum.array() == 0).count() << std::endl;

    // voxelize sanity check
    // S = (S.array() != 0).select(-S, S);
    // igl::copyleft::marching_cubes(S, GV, side(0), side(1), side(2), -1.0, V_voxels_mc, F_voxels_mc);
    // igl::writeSTL("output_marching_cubes.stl", V_voxels_mc, F_voxels_mc);

    std::cout << "views: " << std::endl;
    for(int i = 0; i < views.rows(); i++)
    {
      std::cout << views.row(i) << std::endl;
    }

  glfwDestroyWindow(window);
  glfwTerminate();

  /*
  // transform back meshes
  Eigen::Vector3f m = V.colwise().maxCoeff();
  V *= scale_factor;
  V.rowwise() += translation;

  // std::cout << (V.colwise().maxCoeff()-V.colwise().minCoeff()).maxCoeff() << std::endl;
  V_voxels *= scale_factor;
  V_voxels.rowwise() += translation;

  GV *= scale_factor;
  GV.rowwise() += translation;

  views *= scale_factor;
  views.rowwise() += translation;
  views.rowwise() += Eigen::RowVector3f(0,0,-translation(2)+scale_factor);

  // Create a libigl Viewer object
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V_voxels.cast <double> (), F_voxels);
  viewer.data().set_face_based(true);
  viewer.data().set_colors(Eigen::RowVector3d(146,197,222)/255.);
  // viewer.append_mesh();
  // viewer.data().set_mesh(V.cast <double> (), F_int);
  // viewer.data().set_face_based(true);
  // viewer.data().set_colors(Eigen::RowVector3d(244,165,130)/255.);

  // viewer.data().set_points(GV.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
  // viewer.data().point_size = 10;

  // viewer.data().set_points(views.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
  viewer.launch();
  */  

  return EXIT_SUCCESS;

}