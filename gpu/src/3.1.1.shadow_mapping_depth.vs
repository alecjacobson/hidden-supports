#version 410 core
layout (location = 0) in mediump vec3 aPos;

// const int MAX_VIEWS = 100; 

uniform mediump mat4 light_proj;
uniform mediump mat4 proj;
uniform mediump mat4 model;
// uniform mat4 light_views[MAX_VIEWS];
uniform mediump mat4 light_view;
uniform mediump mat4 view;

void main()
{

    gl_Position = light_proj * light_view * model * vec4(aPos, 1.0);
    //gl_Position = proj * view * model * vec4(aPos, 1.0);
    //gl_Position = proj * model * vec4(aPos, 1.0);
}