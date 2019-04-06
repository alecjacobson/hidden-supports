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
GLuint FBO,FBO_render,FBO_large;
GLuint shadow_map,visibility_map,large_visibility;
// GLuint fbo_render,visible_slice;

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


  std::cout<<R"(
  Usage:
    [Click and drag]  to orbit view
    [Scroll]  to translate view in and out
    A,a  toggle animation
    L,l  toggle wireframe rending
    Z,z  reset view to look along z-axis
  )";

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

  std::cout << "original max and min" << std::endl;
  std::cout << V.colwise().maxCoeff() << std::endl;
  std::cout << V.colwise().minCoeff() << std::endl;
  Eigen::RowVector3f translation = V.colwise().mean();
  // translation = Eigen::RowVector3f(translation(0), translation(1), 0.0);
  V.rowwise() -= translation;
  std::cout << "translation " << translation << std::endl;
  std::cout << "translated max and min" << std::endl;
  std::cout << V.colwise().maxCoeff() << std::endl;
  std::cout << V.colwise().minCoeff() << std::endl;
  float scale_factor = (V.colwise().maxCoeff()-V.colwise().minCoeff()).maxCoeff();
  std::cout << "scale factor " << scale_factor << std::endl;
  V /= scale_factor;
  std::cout << "scaled max and min" << std::endl;
  std::cout << V.colwise().maxCoeff() << std::endl;
  std::cout << V.colwise().minCoeff() << std::endl;

  Eigen::Vector3f centroid;
  igl::centroid(V,F,centroid);
  std::cout << centroid << " CENTROID " << std::endl;
  // look_at(Eigen::Vector3f(0,0,-1), Eigen::Vector3f(0,1,0), centroid, view);

  // // make voxel grid
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
  std::cout << ratio << std::endl;

  max_z = V.col(2).maxCoeff();
  min_z = V.col(2).minCoeff();
  std::cout << max_z << std::endl << min_z << std::endl;
  z_range = std::abs(max_z - min_z);
  int number_of_slices = side(2);
  step = z_range / number_of_slices;

  // glfwSetWindowShouldClose(window,true);
  // window = glfwCreateWindow(w, h, "visibility", NULL, NULL);
  // glfwSetWindowPos(window,0,0);
  // glfwMakeContextCurrent(window);

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
	generate_views(bottom_left, top_right, num_views, views);
	std::cout << "generated viewing rays" << std::endl;
  // std::vector<GLfloat*> light_views((int)num_views);

  // precompute light view matrices
  /*
  for(int v = 0; v < views.rows(); v++)
  {
    Eigen::Vector3f viewpoint = views.row(v);
    Eigen::Vector3f l = centroid - viewpoint;
    l.normalize();
    Eigen::Vector3f n = -(viewpoint - Eigen::Vector3f(0,0,1));
    n.normalize();
    float cos_spin_angle = n.dot(l);
    // if x value is negative, then rotate positive angle
    float spin_angle = (viewpoint(0) < 0) ? acos(cos_spin_angle) : -1*acos(cos_spin_angle);

    Eigen::Affine3f light_view = Eigen::Affine3f::Identity() * 
      Eigen::Translation3f(views.row(v));
    light_view.rotate(
      Eigen::AngleAxisf(
        spin_angle,
        Eigen::Vector3f(0,1,0)));

    light_views[v] = light_view.data();
  }
  */

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
  glfwSetCharModsCallback(
    window,
    [](GLFWwindow* window, unsigned int codepoint, int modifier)
    {
      switch(codepoint)
      {
        case 'O':
        case 'o':
          // Q_V.col(2).array() += 0.1;
          // mesh_to_vao(q_prog_id, Q_V, Q_F, Q_N, Q_TC, Q_VAO);
          z_slice += step;
          std::cout << "quad z " << z_slice << std::endl;
          q_model << ratio, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, z_slice,
                  0, 0, 0, 1;
          break;
        case 'P':
        case 'p':
          // Q_V.col(2).array() -= 0.1;
          // mesh_to_vao(q_prog_id, Q_V, Q_F, Q_N, Q_TC, Q_VAO);
          z_slice -= step;
          std::cout << "quad z " << z_slice << std::endl;
          q_model << ratio, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, z_slice,
                  0, 0, 0, 1;
          break;
        case 'D':
        case 'd':
          view.matrix().block(0,0,3,3).setIdentity();
          break;
        case 'L':
        case 'l':
          wire_frame ^= 1;
          std::cout << "wireframe = " << wire_frame << std::endl;
          break;
        case 'Z':
        case 'z':
          break;
        case 'X':
        case 'x':
          break;
        case 'A':
        case 'a':
          break;
        case 'S':
        case 's':
          break;
        case 'C':
        case 'c':
          break;
        case 'V':
        case 'v':
          break;
        default:
          std::cout<<"Unrecognized key: "<<(unsigned char) codepoint<<std::endl;
          break;
      }
    });

glfwSetMouseButtonCallback(
    window,
    [](GLFWwindow * window, int button, int action, int mods)
    {
      mouse_down = action == GLFW_PRESS;
    });
glfwSetCursorPosCallback(
  window,
  [](GLFWwindow * window, double x, double y)
  {
    static double mouse_last_x = x;
    static double mouse_last_y = y;
    double dx = x-mouse_last_x;
    double dy = y-mouse_last_y;
    if(mouse_down)
    {
      // // Two axis valuator with fixed up
      // float factor = std::abs(light_view.matrix()(2,3));
      // light_view.rotate(
      //   Eigen::AngleAxisf(
      //     dx*factor/float(w),
      //     Eigen::Vector3f(0,1,0)));
      // light_view.rotate(
      //   Eigen::AngleAxisf(
      //     dy*factor/float(h),
      //     light_view.matrix().topLeftCorner(3,3).inverse()*Eigen::Vector3f(1,0,0)));

      // Two axis valuator with fixed up
      // float factor = std::abs(view.matrix()(2,3));
      // view.rotate(
      //   Eigen::AngleAxisf(
      //     dx*factor/float(w),
      //     Eigen::Vector3f(0,1,0)));
      // view.rotate(
      //   Eigen::AngleAxisf(
      //     dy*factor/float(h),
      //     view.matrix().topLeftCorner(3,3).inverse()*Eigen::Vector3f(1,0,0)));
    }
    mouse_last_x = x;
    mouse_last_y = y;
    });
  glfwSetScrollCallback(window,
    [](GLFWwindow * window, double xoffset, double yoffset)
    {
      view.matrix()(2,3) =
        std::min(std::max(view.matrix()(2,3)+(float)yoffset,-100.0f),-0.1f);
      std::cout << "view: " << view.matrix()(2,3) << std::endl;
      
    });

  glEnable(GL_DEPTH_TEST);
  igl::opengl::report_gl_error("enabled\n");

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

  init_shadow_buffer(shadow_map, FBO, q_prog_id, t_w, t_h, "depth");
  igl::opengl::report_gl_error("init shadow buffer\n");

  init_shadow_buffer(visibility_map, FBO_render, render_prog_id, w, h, "color");
  igl::opengl::report_gl_error("init shadow buffer\n");

  init_shadow_buffer(large_visibility, FBO_large, render_prog_id, w*num_views, h, "color");
  igl::opengl::report_gl_error("init large shadow buffer\n");
  std::cout << "SIZE " << w*num_views << " by " << h << std::endl;
 
  // render buffer as color buffer
  // glGenFramebuffers(1,&fbo_render);
  // glBindFramebuffer(GL_FRAMEBUFFER, fbo_render);
  // glGenFramebuffers(1,&fbo_render);
  // glGenRenderbuffers(1, &visible_slice);

  // glBindRenderbuffer(GL_RENDERBUFFER, visible_slice);
  // glRenderbufferStorage(GL_RENDERBUFFER, GL_RED, w, h);
  // glBindFramebuffer(GL_FRAMEBUFFER, fbo_render);

  // attach render buffer to the fbo as color buffer
  // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, visible_slice);
  // bool status = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
  // if(!status)
  //     std::cout << "Could not initialise FBO" << std::endl;
  // else
  //     std::cout << "color FBO ready!" << std::endl;

  float start_time = igl::get_seconds();

  // float factor = std::abs(light_view.matrix()(2,3));
  // light_view.rotate(
  //     Eigen::AngleAxisf(
  //     -0.1*factor/float(w),
  //     Eigen::Vector3f(0,1,0)));

  z_slice = min_z + (number_of_slices/2)*step - step;
  max_z = min_z + (number_of_slices/2)*step;
  Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_values;
  visibility_values.resize(w*h*number_of_slices,1);

  int count = 0;

  // Main display routine
  // while (!glfwWindowShouldClose(window))
  // {
  while(z_slice < max_z)
  {
    std::cout << "Z SLICE NUMBER: " << count << std::endl;
    std::cout << "Z SLICE DEPTH: " << z_slice << std::endl;
    std::cout << "-----" << std::endl;

    igl::opengl::report_gl_error("start of loop\n");

    double tic = igl::get_seconds();

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


      for(int v = 0; v < views.rows(); v++)
      {
          // std::cout << "view " << views.row(v) << std::endl; 
          // look_at(views.row(v), Eigen::Vector3f(0,1,0), centroid, light_view);
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

          // clear screen and set viewport
          glClearColor(1.0,1.0,1.0,0.);
          glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
          igl::opengl::report_gl_error("cleared screen\n");


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

          igl::opengl::report_gl_error("uniforms\n");

          //////////

          bind_map_for_writing(FBO);

          glClear(GL_DEPTH_BUFFER_BIT);

          mesh_to_vao(prog_id, V, F, N, TC,VAO);
          igl::opengl::report_gl_error("bind vao 1\n");
          
          glBindVertexArray(VAO);

          glViewport(0, 0, t_w, t_h);

          // draw elements to texture
          glDrawElements(GL_TRIANGLES, F.size(), GL_UNSIGNED_INT, 0);
          igl::opengl::report_gl_error("draw elements 1\n");

          // Eigen::Matrix< GLubyte,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> shadow;
          // shadow.resize(t_w*t_h,1);

          // glReadPixels(0, 0, t_w, t_h, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, shadow.data());
          // igl::writeDMAT("texture"+std::to_string(v)+".dmat", shadow,true);

          glBindFramebuffer(GL_FRAMEBUFFER, 0);
          igl::opengl::report_gl_error("unbind fbo 0\n");

          // glfwGetFramebufferSize(window, &::w, &::h);
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
          // GLint q_light_source_loc = glGetUniformLocation(q_prog_id,"light_source");
          // glUniform3fv(q_light_source_loc,1,viewpoint.data());
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
          // glUniformMatrix4fv(q_light_view_loc,num_views,GL_FALSE,light_views[0]);
          glUniformMatrix4fv(q_light_view_loc,1,GL_FALSE,light_view.data());

          mesh_to_vao(q_prog_id, Q_V, Q_F, Q_N, Q_TC,Q_VAO);

          igl::opengl::report_gl_error("bind frame buffer for color reading\n");
          bind_map_for_reading(shadow_map, GL_TEXTURE0);
          bind_map_for_writing(FBO_render);

          glBindVertexArray(Q_VAO);
          
          igl::opengl::report_gl_error("bind vao quad\n");
          glDrawElements(GL_TRIANGLES, Q_F.size(), GL_UNSIGNED_INT, 0);
          igl::opengl::report_gl_error("draw elements to offscreen buffer\n");

          Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slice;
                visibility_slice.resize(w*h,1);
          glReadPixels(0, 0, w, h, GL_RED, GL_FLOAT, visibility_slice.data());
          igl::writeDMAT("slice"+std::to_string(v)+".dmat", visibility_slice, true);

          bind_map_for_reading(visibility_map, GL_TEXTURE1);

          GLint x_offset = w*v;
          std::cout << x_offset << std::endl;
          glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, x_offset, 0, w, h);
          igl::opengl::report_gl_error("copy texture\n");

          // if(count < side(2))
          //   visibility_values.block(count*w*h, 0,
          //                           w*h, 1) += visibility_slice;

          // igl::writeDMAT("output_slice_"+std::to_string(count)+".dmat", visibility_slice,true);

          glBindVertexArray(0);

          ///////////

          glFlush();
          glfwSwapBuffers(window);
          igl::opengl::report_gl_error("flush\n");

          {
            glfwPollEvents();
            // In microseconds
            double duration = 1000000.*(igl::get_seconds()-tic);
            const double min_duration = 1000000./60.;
            if(duration<min_duration)
            {
              std::this_thread::sleep_for(std::chrono::microseconds((int)(min_duration-duration)));
            }
          }
          igl::opengl::report_gl_error("poll\n");
         
      }

      Eigen::Matrix< GLfloat,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> visibility_slices;
      visibility_slices.resize(num_views*w*h,1);

      bind_map_for_reading(large_visibility,GL_TEXTURE2);
      glReadPixels(0, 0, w*num_views, h, GL_RED, GL_FLOAT, visibility_slices.data());
      igl::opengl::report_gl_error("read pixels\n");

      igl::writeDMAT("slices.dmat", visibility_slices, true);


      //if(count < side(2))
      //   visibility_values.block(count*w*h, 0,
      //                           w*h, 1) += visibility_slice;

      // glFlush();
      // glfwSwapBuffers(window);

      // {
      //   glfwPollEvents();
      //   // In microseconds
      //   double duration = 1000000.*(igl::get_seconds()-tic);
      //   const double min_duration = 1000000./60.;
      //   if(duration<min_duration)
      //   {
      //     std::this_thread::sleep_for(std::chrono::microseconds((int)(min_duration-duration)));
      //   }
      // }

      z_slice += step;
      count++;
      
    }

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

  double isovalue = 0.8;
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

  glDeleteVertexArrays(1,&VAO);
  glDeleteFramebuffers(1,&FBO);

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