#version 410 core
layout (location = 0) in mediump vec3 aPos;
layout (location = 1) in mediump vec2 aTexCoords;

uniform mediump mat4 model;
uniform mediump mat4 q_model;
uniform mediump mat4 proj;
uniform mediump mat4 light_proj;
uniform mediump mat4 view;
uniform mediump mat4 light_view;

uniform mediump vec3 light_source;

out mediump vec2 tex_coord_out;
out mediump vec4 pos_cam_space;
out mediump vec4 pos_light_space;

void main()
{
    tex_coord_out = aTexCoords;
    vec4 pos = proj * view * q_model * vec4(aPos, 1.0);
    pos_cam_space = pos;
    //quad_pos = 
    pos_light_space = light_proj * light_view * q_model * vec4(aPos, 1.0);
    gl_Position = pos;
}