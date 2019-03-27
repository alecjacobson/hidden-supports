// Input:
//   N  3D unit normal vector
// Outputs:
//   T  3D unit tangent vector
//   B  3D unit tangent vector
void tangent(in vec3 N, out vec3 T, out vec3 B)
{
  vec3 x = vec3(1,0,0);
  T = normalize(x-dot(N,x)*N);
  B = cross(N,T);
}
