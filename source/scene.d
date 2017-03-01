module scene;
import opencl;
import globals;

struct Triangle {
  cl_float3 A, B, C;
  cl_uint material_index;
  cl_uint bufsize1, bufsize2, bufsize3;
}

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
  Triangle[] vertices;
  this ( string _name ) {
    name = _name;
  }
}

struct Camera {
  cl_float3 position, direction;
  cl_int2 dimensions;
}

auto Construct_Camera(float[3] pos, float[3] dir, int[2] dim) {
  cl_int2 dimensions;
  dimensions.x = dim[0]; dimensions.y = dim[1];
  return Camera(
    To_CLFloat3(pos), To_CLFloat3(dir), dimensions
  );
}

auto To_CLFloat3(float[3] a) {
  cl_float3 vec;
  vec.x = a[0]; vec.y = a[1]; vec.z = a[2];
  return vec;
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
  { // --- create model from JSON ---
    static import std.file;
    import functional;
    auto list = std.file.dirEntries(scene_name, std.file.SpanMode.breadth)
                   .filter!(n => n.isFile)
                   .filter!(n => cast(string)(n)[$-4..$] == ".aoq");
    foreach ( fname; list ) {
      auto data = read(fname).to!string.split("\n") // split to array
                  .filter!(n => n.length > 2).array; // remove blank lines

      auto RIndex() {
        auto res = data.filter!(n => n.length > 6)
                       .filter!(n => n[0..6] == "index ").array;
        if ( res.length > 0 ) return res[0][6..$].to!int;
        return 1;
      }
      int index = RIndex();
      writeln("Index: ", index);
      string name = data.filter!(n => n[0..2] == "n ").array[0][2..$];
      writeln("READING DATA FOR: ", name);
      float[3][] vertices;
      size_t[4][] faces;
      data.filter!(n => n[0..2] == "v ")
          .map!(n => n[2..$]
            .split(" ")
            .filter!(n => n.length > 0)
            .map!(n => n.to!float))
          .each!(n => vertices ~= [n.array[0 .. 3]]);
      data.filter!(n => n[0..2] == "f ")
          .map!(n => n[2..$].split(" ")
          .filter!(n => n.chomp.length > 0)
          .map!(n => n.to!size_t - index))
          .each!(n => faces ~= [n.array[0 .. 4]]);
      string[] materials = data.filter!(n => n[0..2] == "m ")
                               .map!(n => n[2..$]).array;
      writeln("NAME: ", name);
      writeln("VERTICES: \n", vertices.length, "\n----");
      writeln("FACES: \n", faces.length, "\n----");
      writeln("MATERIALS: \n", materials, "\n----");
      writeln("-------------");
      foreach ( face; faces ) {
        scene.vertices ~= Triangle(
          To_CLFloat3(vertices[face[0]]), To_CLFloat3(vertices[face[1]]),
          To_CLFloat3(vertices[face[2]]),
          cast(uint)material_indices[materials[face[3]]]
        );
      }
      writeln(scene.vertices);
    }
  }
  return cast(immutable)scene;
}
