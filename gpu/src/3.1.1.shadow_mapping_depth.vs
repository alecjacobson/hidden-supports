#version 410 core
layout (location = 0) in vec3 aPos;

// const int MAX_VIEWS = 100; 

uniform mat4 light_proj;
uniform mat4 proj;
uniform mat4 model;
// uniform mat4 light_views[MAX_VIEWS];
uniform mat4 light_view;
uniform mat4 view;

void main()
{

    gl_Position = light_proj * light_view * model * vec4(aPos, 1.0);
    //gl_Position = proj * view * model * vec4(aPos, 1.0);
    //gl_Position = proj * model * vec4(aPos, 1.0);
}