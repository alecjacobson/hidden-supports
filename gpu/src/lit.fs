// Add (hard code) an orbiting (point or directional) light to the scene. Light
// the scene using the Blinn-Phong Lighting Model.
//
// Uniforms:
uniform mat4 view;
uniform mat4 proj;

uniform float animation_seconds;


// Inputs:
in vec3 n_fs_in;
in vec4 pos_fs_in; 
// Outputs:
out vec3 color;
// expects: blinn_phong, identity
void main()
{
  vec3 blue = vec3(0.2,0.3,0.8);
  vec3 gray = vec3(0.5,0.45,0.5);
  
  float light_freq = 0.125*2.0*M_PI;
  vec3 l = mat3(view)*normalize(vec3(cos(light_freq*animation_seconds),0.8,(sin(light_freq*animation_seconds))));
  //vec3 l = mat3(view)*normalize(vec3(15,15,15));
  vec3 n = normalize(n_fs_in);
  mat4 M = identity();
  vec4 view_pos_fs_in = proj*view*M*pos_fs_in;
  vec3 v = normalize(-view_pos_fs_in.xyz);
  color = blinn_phong(
    blue,
    gray,
    vec3(1,1,1),
    100.0,
    n,
    v,
    l);
}
