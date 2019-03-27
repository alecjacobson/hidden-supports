// Set the pixel color to blue or gray depending on is_moon.
//
// Uniforms:
uniform bool is_moon;
// Outputs:
out vec3 color;
void main()
{
  vec3 blue = vec3(0.2,0.3,0.8);
  vec3 gray = vec3(0.5,0.45,0.5);
  color = mix(blue,gray,float(is_moon));
}
