#version 410 core
layout (location = 0) in highp vec3 aPos;
layout (location = 1) in highp vec2 aTexCoords;

uniform highp mat4 model;
uniform highp mat4 q_model;
uniform highp mat4 proj;
uniform highp mat4 light_proj;
uniform highp mat4 view;
uniform highp mat4 light_view;

uniform highp vec3 light_source;

out highp vec2 tex_coord_out;
out highp vec4 pos_cam_space;
out highp vec4 pos_light_space;

void main()
{
    tex_coord_out = aTexCoords;
    vec4 pos = proj * view * q_model * vec4(aPos, 1.0);
    pos_cam_space = pos;
    //quad_pos = 
    pos_light_space = light_proj * light_view * q_model * vec4(aPos, 1.0);
    gl_Position = pos;
}