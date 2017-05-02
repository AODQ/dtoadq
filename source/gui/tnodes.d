module gui.tnodes;
import derelict.imgui.imgui : ImVec2;
import globals, gui.gui, std.variant : Variant;

// ----------------------------------------------------------------------------
// --------------------------- public access ----------------------------------
auto RNodes ( ) { return g_nodes; }
auto RNode ( int ID ) in {
  assert(ID in g_nodes);
} body {
  return g_nodes[ID];
}

auto RNodeConnection ( int ID   ) { return g_node_connections[ID]; }
auto RNodeType       ( string name ) in {
  assert(name in g_node_types, "RNodeType from invalid name: " ~ name);
} body { return g_node_types[name];     }

auto New_Node ( string name, ImVec2 vec ) {
  writeln("Adding node!");
  g_nodes[g_node_id_counter] = new Node(name.RNodeType, g_node_id_counter, vec);
  writeln(g_nodes);
  ++ g_node_id_counter;
}
auto New_Node ( string name, float x, float y ) {
  New_Node(name, ImVec2(x, y));
}

auto Add (inout ImVec2 x, inout ImVec2 y) {return ImVec2(x.x + y.x, x.y + y.y);}
auto Sub (inout ImVec2 x, inout ImVec2 y) {return ImVec2(x.x - y.x, x.y - y.y);}

NodeConnection RHovering_Node ( ImVec2 offset, out ImVec2 origin ) {
  foreach ( ref node; g_nodes ) {
    ImVec2 node_origin = Add(node.origin, offset);
    foreach ( it, ref con; node.RConnection_Descs ) {
      if ( !con.RHovering(offset) ) continue;
      origin = Add(node_origin, con.ROrigin);
      return new NodeConnection(node.RID, cast(int)it);
    }
  }
  return null;
}

unittest {
  New_Node("Add", 0.0, 0.0);
  assert(g_nodes.length > 0 && g_nodes[0] !is null);
  auto node = g_nodes[0];
  assert(node.RID == 0);
  assert(node.name    == "Add");
  assert(node.subnode_data.type.name == "Add");
}
// ----------------------------------------------------------------------------
// --------------------------- globals ----------------------------------------
immutable float Node_slot_radius = 4.0f;
immutable ImVec2 Node_window_padding = ImVec2(8.0f, 8.0f);
private:
int                 g_node_id_counter, g_node_connections_counter;
Node          [int] g_nodes;
NodeConnection[int] g_node_connections;
NodeType[string] g_node_types;

static this() {
  g_node_types = [
  "Add" : new NodeType("Add", [SubnodeType("In0", SubnodeDescription.Varying),
                               SubnodeType("In1", SubnodeDescription.Varying)],
                              [SubnodeType("Out", SubnodeDescription.Varying)]),
  ];
}

// ----------------------------------------------------------------------------
// --------------------------- structures -------------------------------------
struct ConnectorBase(T) {
  T[] inputs, outputs;
}

enum ConnectionType { Input, Output };
public class NodeConnection {
  int in_node_id, out_node_id;
  int in_subnode_id, out_subnode_id;

  this ( int in_node_id_, int in_subnode_id_ ) {
    in_node_id = in_node_id_;
    in_subnode_id = in_subnode_id_;
  }
}

public class Node {
private:
  ImVec2 origin, size;
  int id;
  string name;
  SubnodeData subnode_data;
public:
  this ( NodeType node_type, int id_, ImVec2 origin_ ) {
    id = id_;
    subnode_data = new SubnodeData(node_type);
    name = node_type.name;
    origin = origin_;
    // -- get size --
    import std.algorithm : max;
    size_t length_in, length_out;
    foreach ( i; subnode_data.data.inputs )
      length_in = max(length_in, i.name.length);
    foreach ( i; subnode_data.data.outputs )
      length_out = max(length_out, i.name.length);
    writeln("length in: ", length_in);
    writeln("length out: ", length_out);
    writeln(subnode_data.data.outputs);
    size = ImVec2((length_in + length_out)*20.0f + 30.0f,
                  max(subnode_data.data.inputs.length,
                      subnode_data.data.outputs.length)*30.0f + 50.0f);
  }

  auto RID   ( ) inout { return id; }
  auto RName ( ) inout { return name; }
  auto ROrigin ( ) inout { return origin; }
  auto RSize ( ) inout { return size; }

  auto SOrigin ( ImVec2 origin_ ) {
    import std.math : fmax;
    origin = origin_;
  }

  auto RInput_Connection_Descs ( ) inout {
    return subnode_data.type.subnode_descs.inputs;
  }
  auto ROutput_Connection_Descs ( ) inout {
    return subnode_data.type.subnode_descs.outputs;
  }
  auto RConnection_Descs ( ) inout {
    return RInput_Connection_Descs ~ ROutput_Connection_Descs;
  }
  auto RSubnode_Data ( ) inout { return subnode_data; }
}

public enum SubnodeDataType { Data, Connection }

class SubnodeData {
  public:
  struct SubnodeDataContainer {
    union DataUnion {
      Variant data;
      int connection_id;
    }
    SubnodeDataType type;
    string name;
    DataUnion data;
    // @disable this;
    this ( string name_, Variant data_ ) {
      name = name_; type = SubnodeDataType.Data;
      data.data = data_;
      writeln("RESNAME: ", name);
    }
    this ( string name_, int connection_id_ ) {
      name = name_; type = SubnodeDataType.Connection;
      data.connection_id = connection_id_;
    }
  }

  NodeType type;
  ConnectorBase!SubnodeDataContainer data;

  static auto Create_Default_Data ( SubnodeDescription info ) {
    final switch ( info ) with ( SubnodeDescription ) {
      case Float3: case Colour: return Variant(gln.vec3(0.0f));
      case Float2:              return Variant(gln.vec2(0.0f));
      case Float:               return Variant(0.0f);
      case Int: case Varying:   return Variant(cast(int)0);
      case String:              return Variant("");
      case Texture:             return Variant("N/A");
      case BRDF:                return Variant("N/A");
    }
  }
public:
  this ( NodeType type_ ) {
    import functional;
    type = type_;
    data.inputs = type_.subnode_descs.inputs.map!((SubnodeType n) {
      return SubnodeDataContainer(n.name, Create_Default_Data(n.description));
    }).array;
    data.outputs = type_.subnode_descs.outputs.map!((SubnodeType n) {
      return SubnodeDataContainer(n.name, Create_Default_Data(n.description));
    }).array;
    data.inputs.length  = type_.RInput_Length;
    data.outputs.length = type_.ROutput_Length;
  }
}

enum SubnodeDescription {
  Colour, Float3, Float2, Float, Int, Varying, String,

  Texture, BRDF
};

struct SubnodeType {
  string             name;
  SubnodeDescription description;
  auto ROrigin ( ) {
    return ImVec2(name.length*13, 13);
  }

  bool RHovering ( ImVec2 offset ) {
    auto p = Sub(gdRMousePos, Add(offset, ROrigin));
    return (p.x*p.x + p.y*p.y) < (Node_slot_radius*Node_slot_radius);
  }
}

class NodeType {
private:
  string name;
  ConnectorBase!SubnodeType subnode_descs;
public:
  this ( string name_, SubnodeType[] inputs, SubnodeType[] outputs ) {
    name = name_;
    subnode_descs = ConnectorBase!SubnodeType(inputs, outputs);
  }

  auto RInput_Length  ( ) inout { return subnode_descs.inputs.length;  }
  auto ROutput_Length ( ) inout { return subnode_descs.outputs.length; }
}


