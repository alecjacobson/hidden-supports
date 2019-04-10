#version 410 core
out mediump vec4 FragColor;

in mediump vec2 tex_coord_out;
in mediump vec4 pos_cam_space;
in mediump vec4 pos_light_space;

uniform mediump sampler2D shadow_map;
uniform mediump sampler2D visibility_map;

uniform mediump float near_plane;
uniform mediump float far_plane;

uniform mediump int index;

uniform mediump mat4 proj;
uniform mediump mat4 light_proj;
uniform mediump mat4 view;
uniform mediump mat4 light_view;

uniform mediump vec3 light_source;



float epsilon = 1e-3; 

bool essentiallyEqual(float a, float b)
{
    return abs(a - b) <= ( (abs(a) > abs(b) ? abs(b) : abs(a)) * epsilon);
}

// required when using a perspective projection matrix
float LinearizeDepth(float depth)
{
    float z = depth * 2.0 - 1.0; // Back to NDC 
    return (2.0 * near_plane * far_plane) / (far_plane + near_plane - z * (far_plane - near_plane));	
}

float ShadowCalculation(vec4 light_space_position, vec4 cam_space_position)
{

    // project point from camera space into light space
    vec3 projected_quad = cam_space_position.xyz / cam_space_position.w;
    projected_quad = projected_quad * 0.5 + 0.5;

    // perform perspective divide
    vec4 light_source_light_space = light_proj * light_view * vec4(light_source,1.0);
    vec3 projected_light = light_source_light_space.xyz / light_source_light_space.w;
    // projected_light = projected_light * 0.5 + 0.5;
    vec3 projected_coords = light_space_position.xyz / light_space_position.w;

    float coords_to_light = distance(light_space_position.xyz, projected_light);

    // transform from [-1,1] to [0,1] range
    projected_coords = projected_coords * 0.5 + 0.5;
    // get closest depth value from light's perspective 
    // (using [0,1] range light_space_position as coords)
    float closest = texture(shadow_map, projected_coords.xy).r;
    //closest = LinearizeDepth(closest) ;// / far_plane;

    // calculate whether object is in shadow
    // if in shadow, then FILL the voxel
    // (distance from projected coordinate to light source) > closest
    float shadow = (projected_coords.z > closest)
                    /*|| essentiallyEqual(projected_quad.z,closest) 
                    || (projected_quad.z + 0.35 < closest) */
                        ? 1.0 : 0.0;

    // keep the shadow at 0.0 when outside the far_plane region of the light's frustum.
    if(projected_coords.z > 1.0)
        shadow = 0.0;
        
    // return texture(visibility_map, projected_coords.xy).r;
    return shadow;

}


void main()
{
    vec3 color = vec3(1,1,1);
    float in_shadow = ShadowCalculation(pos_light_space, pos_cam_space);                      
    vec3 lighting = (in_shadow) * color;

    vec3 projected_coords = pos_light_space.xyz / pos_light_space.w;
    projected_coords = projected_coords * 0.5 + 0.5;

    if (index != 0)
    {
        FragColor = vec4(lighting, 1.0) + texture(visibility_map, tex_coord_out);
    }
    else
    {
        FragColor = vec4(lighting, 1.0);
    }
}