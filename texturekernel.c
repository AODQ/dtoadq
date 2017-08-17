__constant float PI = 3.1415926535f;

//---MAP GEOMETRY INSERTION POINT---
//%MAPFUNCDECLARATIONS
//----------------------------------
//%MAPFUNCDEFINITIONS
//----------------------------------

float4 TextureMap ( float2 origin ) {
  //---MAP INSERTION POINT---
  //%MAPINSERT
  //-------------------------
}

__kernel void Texture_kernel (
                              __global float* img_buffer, __global float* dimensions
                              ) {
  float dims = *dimensions;
  int2 out = (int2)(get_global_id(0), get_global_id(1));
  float4 results = TextureMap((float2)(out.x/dims, out.y/dims));
  int index = dims*out.y*4 + out.x*4;
  img_buffer[index+0] = results.x;
  img_buffer[index+1] = results.y;
  img_buffer[index+2] = results.z;
  img_buffer[index+3] = results.w;
}
