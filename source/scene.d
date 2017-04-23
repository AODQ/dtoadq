module scene;
import opencl;
import globals;

struct Material {
  cl_float3 base_colour;
  cl_float odd_buffer;
  cl_float metallic, subsurface, specular, roughness, specular_tint,
           anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
           emission;
}


struct Camera {
  cl_float3 position, lookat, up;
  cl_int2 dimensions;
  cl_float fov;
}

auto Construct_Camera(float[3] pos, float[3] dir, int[2] dim) {
  cl_int2 dimensions;
  dimensions.x = dim[0]; dimensions.y = dim[1];
  return Camera(
    To_CLFloat3(pos), To_CLFloat3(dir), To_CLFloat3([0.0f, 1.0f, 0.0]),
    dimensions,
    60.0f
  );
}

void Print_Material(ref inout(Material) m) {
  writeln("Base Colour: ",     m.base_colour     );
  writeln("Metallic   : ",     m.metallic        );
  writeln("Subsurface: ",      m.subsurface      );
  writeln("Specular: ",        m.specular        );
  writeln("Roughness: ",       m.roughness       );
  writeln("Specular_tint: ",   m.specular_tint   );
  writeln("Anisotropic: ",     m.anisotropic     );
  writeln("Sheen: ",           m.sheen           );
  writeln("Sheen_tint: ",      m.sheen_tint      );
  writeln("Clearcoat: ",       m.clearcoat       );
  writeln("Clearcoat_gloss: ", m.clearcoat_gloss );
  writeln("Emission: ",        m.emission        );
}
