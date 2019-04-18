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
#include <igl/remove_unreferenced.h>
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
#include "write_pgm.h"

// int w=512,h=301;
int w=100,h=90;
int d=100;
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
float light_top = tan((40./2.)*M_PI/180.)*near;
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
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

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
  z_range = std::abs(max_z - min_z);
  int number_of_slices = side(2);
  step = z_range / number_of_slices;

  // Projection and modelview matrices
  perspective(-light_right, light_right, -light_top, light_top, near, far, light_proj);
  // Orthographic projection matrix
  orthographic(-right, right, -top, top, near, far, proj);

  // get view rays
  float num_views = 500.0;
  // Eigen::Vector3f bottom_left = V.colwise().minCoeff();
  Eigen::Vector3f bottom_left(-0.5,0,-1);  
  // Eigen::Vector3f top_right = V.colwise().maxCoeff();
  Eigen::Vector3f top_right(0.5,0.4,-1);
  Eigen::MatrixXf views;
  std::vector<float> z_vals{-1, 1};
  generate_views(bottom_left, top_right, num_views, z_vals, views);
  // views.resize(2,3);
  // views.row(0) = Eigen::Vector3f(0,0.4,-1);
  // views.row(1) = Eigen::Vector3f(0,0.4,1);
	std::cout << "generated viewing rays" << std::endl;
  std::cout << "views: " << std::endl;
  for(int i = 0; i < views.rows(); i++)
  {
    std::cout << views.row(i) << std::endl;
  }
  igl::writeDMAT("../results/views.dmat", views, true);

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
      Eigen::MatrixXf V_voxels;
      Eigen::MatrixXi F_voxels;

  while(z_slice < max_z-step)
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

      init_shadow_buffer(shadow_map, FBO, GL_TEXTURE0, t_w, t_h, "depth");
      igl::opengl::report_gl_error("init shadow buffer\n");

      init_shadow_buffer(visibility_map_odd, FBO_render_odd, GL_TEXTURE1, w, h, "color");
      igl::opengl::report_gl_error("init shadow buffer\n");

      init_shadow_buffer(visibility_map_even, FBO_render_even, GL_TEXTURE1, w, h, "color");
      igl::opengl::report_gl_error("init shadow buffer\n");

      for(int v = 0; v < views.rows(); v++)
      {
        

        // std::cout << "view tic " << v << " " << tictoc() << std::endl;

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

    glDeleteVertexArrays(1,&VAO);
    glDeleteFramebuffers(1,&FBO_large);

    std::cout << "size of vector: " << visibility_values.rows() << std::endl;
    Eigen::MatrixXi vv = visibility_values.cast <int> ();
    Eigen::VectorXi S(Eigen::Map<Eigen::VectorXi>(vv.data(), vv.cols()*vv.rows()));      

    igl::writeDMAT("../results/visibilities.dmat", S, true);

    Eigen::MatrixXf V_voxels_mc;
    Eigen::MatrixXi F_voxels_mc;
    std::cout << side << std::endl;

    for(int isovalue = 1; isovalue <= views.rows(); isovalue+=10)
    {
      std::cout << "isovalue " << isovalue << std::endl;
      std::stringstream ss;
      ss << std::setfill('0') << std::setw(3) << isovalue;

      // Eigen::VectorXf S_iso = (S.array() >= isovalue).select(1, S.array()-S.array());
      // igl::writeDMAT("visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

      // voxelize
      Eigen::MatrixXf V_hex;
      Eigen::MatrixXi F_hex;
      Eigen::MatrixXf V_quad;
      Eigen::MatrixXi F_quad;
      Eigen::MatrixXi I_quad;

      make_voxels_from_visibility(S, GV, side, (double)isovalue, V_voxels, F_voxels);

      // make_hex_from_visibility(S, GV, side, isovalue, V_hex, F_hex);

      // igl::writeDMAT("V_hex_"+ss.str()+".dmat", V_hex, true);
      // igl::writeDMAT("F_hex_"+ss.str()+".dmat", F_hex, true);
      
      // extract_voxel_surface(V_hex, F_hex, V_quad, F_quad, I_quad);
     
      Eigen::VectorXi S_iso = (S.array() >= isovalue).select(isovalue, S.array()-S.array());
      igl::writeDMAT("../results/visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

      // igl::writeOBJ("../results/output_voxels_"+ss.str()+".obj", V_hex, F_hex);
      // igl::writeOBJ("output_voxels_"+ss.str()+".obj", V_quad, F_quad);
      write_pgm("output_vol_"+ss.str()+".pgm3d", S,side,isovalue);

    }

  glfwDestroyWindow(window);
  glfwTerminate();

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
  // views.rowwise() += Eigen::RowVector3f(0,0,-translation(2)+scale_factor);

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

  viewer.data().set_points(views.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
  viewer.data().point_size = 10;
  viewer.launch();

  return EXIT_SUCCESS;

}
