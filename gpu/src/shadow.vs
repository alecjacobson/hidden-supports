#version 410 core
layout (location = 0) in vec3 aPos;

out vec4 light_pos;

uniform mat4 proj;
uniform mat4 light_proj;
uniform mat4 model;
uniform mat4 view;
uniform mat4 light_view;

void main()
{
    light_pos = light_proj * light_view * model * vec4(aPos, 1.0);
    gl_Position = proj * view * model * vec4(aPos, 1.0);
}