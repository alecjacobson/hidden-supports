#version 410 core

in vec4 light_pos;

out vec4 FragColor;

uniform sampler2D shadow_map;
uniform float near_plane;
uniform float far_plane;

float ShadowCalculation(vec4 light_space_position)
{
    // perform perspective divide
    vec3 projected_coords = light_space_position.xyz / light_space_position.w;
    // transform to [0,1] range
    projected_coords = projected_coords * 0.5 + 0.5;
    // get closest depth value from light's perspective 
    // (using [0,1] range light_space_position as coords)
    float closest = texture(shadow_map, projected_coords.xy).r; 
    // get depth of current fragment from light's perspective
    float current = projected_coords.z;

    // check whether current frag pos is in shadow
    float shadow = current > closest  ? 1.0 : 0.0;
    // PCF

    // float shadow = 0.0;
    // vec2 texel_size = 1.0 / textureSize(shadow_map, 0);
    // for(int x = -1; x <= 1; ++x)
    // {
    //     for(int y = -1; y <= 1; ++y)
    //     {
    //         float pcf_depth = texture(shadow_map, projected_coords.xy + vec2(x, y) * texel_size).r; 
    //         shadow += current > pcf_depth  ? 1.0 : 0.0;        
    //     }    
    // }
    // shadow /= 9.0;
    
    // keep the shadow at 0.0 when outside the far_plane region of the light's frustum.
    // if(projected_coords.z > 1.0)
    //     shadow = 0.0;
        
    return shadow;
}


// required when using a perspective projection matrix
float LinearizeDepth(float depth)
{
    float n = near_plane;// * 100;
    float z = depth * 2.0 - 1.0; // Back to NDC 
    return (2.0 * n * far_plane) / (far_plane + n - z * (far_plane - n));	
}

void main()
{
    vec3 color = vec3(1,1,1);
    // calculate shadow
    float shadow = ShadowCalculation(light_pos);        
    vec3 lighting = (1.0 - shadow) * color;
    
    FragColor = vec4(lighting, 1.0);
}