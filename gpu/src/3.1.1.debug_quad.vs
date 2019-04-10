#version 410 core
layout (location = 0) in  vec3 aPos;
layout (location = 1) in  vec2 aTexCoords;

uniform  mat4 model;
uniform  mat4 q_model;
uniform  mat4 proj;
uniform  mat4 light_proj;
uniform  mat4 view;
uniform  mat4 light_view;

uniform  vec3 light_source;

out  vec2 tex_coord_out;
out  vec4 pos_cam_space;
out  vec4 pos_light_space;

void main()
{
    tex_coord_out = aTexCoords;
    vec4 pos = proj * view * q_model * vec4(aPos, 1.0);
    pos_cam_space = pos;
    //quad_pos = 
    pos_light_space = light_proj * light_view * q_model * vec4(aPos, 1.0);
    gl_Position = pos;
}