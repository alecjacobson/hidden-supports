// Set the pixel color using Blinn-Phong shading (e.g., with constant blue and
// gray material color) with a bumpy texture.
// 
// Uniforms:
uniform mat4 view;
uniform mat4 proj;
uniform float animation_seconds;
uniform bool is_moon;
// Inputs:
//                     linearly interpolated from tessellation evaluation shader
//                     output
in vec3 sphere_fs_in;
in vec3 normal_fs_in;
in vec4 pos_fs_in; 
in vec4 view_pos_fs_in; 
// Outputs:
//               rgb color of this pixel
out vec3 color;
// expects: model, blinn_phong, bump_height, bump_position,
// improved_perlin_noise, tangent
void main()
{
  vec3 blue = vec3(0.2,0.3,0.8);
  vec3 gray = vec3(0.5,0.45,0.5);
  // Rotate light around the scene at a frequency of 1 revolutions per 8 seconds
  float light_freq = 0.125*2.0*M_PI;
  vec3 l = mat3(view)*normalize(vec3(cos(light_freq*animation_seconds),0.8,(sin(light_freq*animation_seconds))));
  vec3 v = normalize(-view_pos_fs_in.xyz);

  vec3 s = normalize( sphere_fs_in );
  vec3 p = normalize( sphere_fs_in );
  vec3 n = normalize( sphere_fs_in );
  float b = bump_height( is_moon , s );
  vec3 bp = bump_position(is_moon, s );

  float eps = 0.0001;
  vec3 T;
  vec3 B;
  tangent(n,T,B);
  vec3 pT = p+eps*T;
  vec3 pB = p+eps*B;
  n = normalize(cross(
    bump_position(is_moon, pT)-bp,
    bump_position(is_moon, pB)-bp));
  mat4 M = model(is_moon,animation_seconds);
  n = normalize((view*M*vec4(n,0)).xyz);

  color = mix(blue,gray,float(is_moon));
  color = vec3(0.5,0.5,0.5) + n*vec3(0.5,0.5,0.5);
  color = blinn_phong(
    clamp(1+5*b,0,1)*mix(blue,gray,float(is_moon)),
    clamp(1+5*b,0,1)*mix(blue,gray,float(is_moon)),
    vec3(1,1,1),
    1000.0,
    n,
    v,
    l);
}
