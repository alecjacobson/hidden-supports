#version 330 core
//in vec3 n_fs_in;
out vec4 color;

in vec2 TexCoordOut;
uniform sampler2D gShadowMap;

const highp vec4 packFactors = vec4(1.0, 255.0, 65025.0, 16581375.0);
const highp vec4 cutoffMask  = vec4(1.0/255.0,1.0/255.0,1.0/255.0,0.0);

void main()
{
    highp vec4 packedVal = vec4(fract(packFactors * gl_FragCoord.z));
    float ndcDepth = 
      (2.0 * gl_FragCoord.z - gl_DepthRange.near - gl_DepthRange.far) /
      (gl_DepthRange.far - gl_DepthRange.near);
    float clipDepth = ndcDepth / gl_FragCoord.w;
    //color = vec4((clipDepth * 0.5) + 0.5); 
    //color = packedVal - packedVal.yzww * cutoffMask;

    float depth = texture(gShadowMap, TexCoordOut).x;
    //depth = 1.0 - (1.0 - depth);// * 10.0;
    color = vec4(vec3(depth), 1.0);
    //color = vec4(0.5,0.5,0.5,1.0) * texture(gShadowMap, TexCoordOut);
    //color = vec4(gl_FragCoord.z);
}