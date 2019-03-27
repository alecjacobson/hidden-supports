#version 330 core
uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;

layout (location = 0) in vec3 position;
//layout (location = 1) in vec3 n_vs_in;
layout (location = 1) in vec2 TexCoord;

//in vec3 position;
//in vec3 n_vs_in;

//out vec3 n_fs_in;
out vec2 TexCoordOut;

void main()
{
  //n_fs_in = n_vs_in;
  //gl_Position = proj * view * model * vec4(position,1.);
  gl_Position = proj * model * vec4(position,1.);
  TexCoordOut = TexCoord;
}
