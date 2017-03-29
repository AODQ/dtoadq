module scene;
import opencl;
import globals;
import cloctree;

struct Material {
  cl_float3 base_colour;
  cl_float odd_buffer;
  cl_float metallic, subsurface, specular, roughness, specular_tint,
           anisotropic, sheen, sheen_tint, clearcoat, clearcoat_gloss,
           emission;
}

class Scene {
public:
  string name;
  Material[] materials;
  OctreeData data;
  this ( string _name ) {
    name = _name;
  }
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


immutable(Scene) Create_Scene(string scene_name) {
  import std.json;
  Scene scene = new Scene(scene_name);
  scene_name = "./scene/"~scene_name~"/";
  size_t[string] material_indices;
  import std.file : read;
  { // --- create materials from JSON ---
    auto contents = read(scene_name~"materials.txt").to!(char[]);
    auto material_json = parseJSON(contents);
    foreach ( material; material_json["materials"].array ) {
      scene.materials ~= Material(
        To_CLFloat3(material["base_colour"].str.to!(float[3])),
        material["metallic"        ].floating,
        material["subsurface"      ].floating,
        material["specular"        ].floating,
        material["roughness"       ].floating,
        material["specular_tint"   ].floating,
        material["anisotropic"     ].floating,
        material["sheen"           ].floating,
        material["sheen_tint"      ].floating,
        material["clearcoat"       ].floating,
        material["clearcoat_gloss" ].floating,
        material["emission"        ].floating
      );
      material_indices[material["name"].str] = scene.materials.length-1;
      writeln("MATERIAL: ", material["name"].str);
      Print_Material(scene.materials[$-1]);
      writeln("----------");
    }
  }
  CLVoxel[] voxels;
  { // --- create model from JSON ---
    static import std.file;
    import functional;
    auto list = std.file.dirEntries(scene_name, std.file.SpanMode.breadth)
                   .filter!(n => n.isFile)
                   .filter!(n => cast(string)(n)[$-4..$] == ".aoq");
    foreach ( fname; list ) {
      auto data = read(fname).to!string.split("\n") // split to array
                  .filter!(n => n.length > 2).array; // remove blank lines

      voxels = data
          .filter!(n => n[0 .. 2] == "v ") // grab voxels
          .map!(n => n[2 .. $] // remove the beginning v
            .split(" ") // split into array of num strings
            .filter!(n => n.length > 0) // remove balanks
            .map!(n => n.to!float).array) // make array of floats
          .map!(n => New_CLVoxel(n[0..3])).array; // make array of voxels
      string[] materials = data.filter!(n => n[0..2] == "m ")
                               .map!(n => n[2..$]).array;
      writeln("NAME: ", fname);
      writeln("VOXELS: \n", voxels.length, "\n----");
      writeln("MATERIALS: \n", materials, "\n----");
    }
  }
  return cast(immutable)scene;
}
