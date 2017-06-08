module gui.tnodes;
import derelict.imgui.imgui : ImVec2, igGetMousePos;
import globals, std.variant : Variant;

// ----------------------------------------------------------------------------
// --------------------------- public access ----------------------------------
auto RNodes ( ) { return g_nodes; }
auto RNode ( int ID ) in {
  assert(ID in g_nodes);
} body {
  return g_nodes[ID];
}

auto Remove_Node ( int ID ) in {
  assert(ID in g_nodes);
} body {
  Clear_Connections(ID);
  writeln("REMOVING: ", ID);
  g_nodes.remove(ID);
  writeln("AFTER: ", g_nodes);
}

auto Clear_Connections ( int ID ) in {
  assert(ID in g_nodes);
} body {
  auto node = g_nodes[ID];
  // we remove both parts of connection, even though we only need to remove one
  // part as this node is deleted
  foreach ( con_tuple; node.RSubnode_Connections ) {
    int con_id = (*con_tuple[1]).data;
    if ( con_id == -1 ) continue;
    auto con = con_id.RNodeConnection;
    con. in_node_id.RNode.Subnode_Connection(con. in_subnode_id, -1, true );
    con.out_node_id.RNode.Subnode_Connection(con.out_subnode_id, -1, false);
    // finally we delete this connection
    g_node_connections.remove(con_id);
  }
}

auto SNode_ID_Counter          ( int id ) { g_node_id_counter          = id; }
auto SNode_Connections_Counter ( int id ) { g_node_connections_counter = id; }
auto RNode_ID_Counter ( ) { return g_node_id_counter; }
auto RNode_Connections_Counter ( ) { return g_node_connections_counter; }
void Clear_Nodes ( ) { g_nodes.clear; }
void Clear_Node_Connections ( ) { g_node_connections.clear; }

auto RSubnode ( int node_id, int subnode_id, bool is_input ) in {
  assert(subnode_id >= 0);
  assert(subnode_id < is_input? RNode(node_id).RInput_Connection_Descs .length
                              : RNode(node_id).ROutput_Connection_Descs.length);
} body {
  auto node = RNode(node_id);
  return is_input ? node.RInput_Connection_Descs[subnode_id] :
                    node.RInput_Connection_Descs[subnode_id];
}

auto RNodeConnection ( int ID   ) { return g_node_connections[ID]; }
auto RNode_Connections( ) { return g_node_connections; }
auto RNodeType       ( string name ) in {
  assert(name in g_node_types, "RNodeType from invalid name: " ~ name);
} body { return g_node_types[name];     }

auto To_Proper_Type ( int type ) {
  switch ( type ) {
    default: assert(0);
    case 0: return SubnodeDescription.Float;
    case 1: return SubnodeDescription.Float2;
    case 2: return SubnodeDescription.Float3;
    case 3: return SubnodeDescription.Int;
  }
}

auto New_Node_From_File ( string filename, ImVec2 vec ) {
  import functional, parser : Parse_File, Create_Model_Info;
  auto info = filename.Parse_File.Create_Model_Info;
  NodeType nodetype = new NodeType(
    info.name,
    info.params.map!(n => SubnodeType(n.label, To_Proper_Type(n.type)))
        .array,
    [SubnodeType("out", SubnodeDescription.Float)],
    ParseNodeType.Function
  );
  auto node = new Node(nodetype, g_node_id_counter, vec);
  g_nodes[g_node_id_counter ++] = node;
  return node;
}

auto New_Node ( string name, ImVec2 vec, int id ) {
  auto node = new Node(name.RNodeType, id, vec);
  g_nodes[id] = node;
  return node;
}
auto New_Node ( string name, ImVec2 vec ) {
  auto node = new Node(name.RNodeType, g_node_id_counter, vec);
  g_nodes[g_node_id_counter ++] = node;
  return node;
}
auto New_Node ( string name, float x, float y ) {
  New_Node(name, ImVec2(x, y));
}

auto Add_Connection ( NodeConnection con, int id ) {
  g_node_connections[id] = new NodeConnection(con);
  // -- set input/output connections --
  g_nodes[con. in_node_id].Subnode_Connection(con. in_subnode_id, con.id, 1);
  g_nodes[con.out_node_id].Subnode_Connection(con.out_subnode_id, con.id, 0);
  return con;
}
auto Add_Connection ( NodeConnection con ) {
  // -- set connection id --
  con.id = g_node_connections_counter ++;
  return Add_Connection(con, con.id);
}

auto Add (inout ImVec2 x, inout ImVec2 y) {return ImVec2(x.x + y.x, x.y + y.y);}
auto Sub (inout ImVec2 x, inout ImVec2 y) {return ImVec2(x.x - y.x, x.y - y.y);}

NodeConnection RHovering_Node ( ImVec2 offset, out ImVec2 origin ) {
  writeln("Offset: ", offset);
  foreach ( ref node; g_nodes ) {
    ImVec2 node_origin = Add(node.origin, offset);
    foreach ( it, ref connection; node.RConnection_Descs ) {
      auto snode = connection[1];
      if ( !snode.RHovering(Sub(node.ROrigin, offset)) ) continue;
      origin = Add(node_origin, snode.ROrigin);
      return new NodeConnection(node.RID, cast(int)it, cast(bool)connection[0]);
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
  import gui.opencl_funcs;
  g_node_types = RNode_Type_Map();
}

// ----------------------------------------------------------------------------
// --------------------------- structures -------------------------------------
struct ConnectorBase(T) {
  T[] inputs, outputs;
}

public auto Text_size = 7.0f;

enum ConnectionType { Input, Output };
public class NodeConnection {
  int in_node_id,    out_node_id,
      in_subnode_id, out_subnode_id;
  int id;

  void Set_In  ( int in_node_id_,  int in_subnode_id_ ) {
    in_node_id   = in_node_id_;  in_subnode_id  = in_subnode_id_;
  }

  void Set_Out ( int out_node_id_, int out_subnode_id_ ) {
    out_node_id  = out_node_id_; out_subnode_id = out_subnode_id_;
  }

  void Set ( int node_id_, int subnode_id_, bool is_input ) {
    (is_input ? &Set_In : &Set_Out)(node_id_, subnode_id_);
  }

  int RNode_ID (bool is_input) { return is_input ? in_node_id : out_node_id; }
  int RSubnode_ID ( bool is_input ) {
    return is_input ? in_subnode_id : out_subnode_id;
  }

  this ( int node_id_, int subnode_id_, bool is_input ) {
    Set(node_id_, subnode_id_, is_input);
  }

  this ( int id_, int in_node_id_,    int out_node_id_,
                  int in_subnode_id_, int out_subnode_id_ ) {
    id = id_;
    in_node_id  = in_node_id_;  in_subnode_id  = in_subnode_id_;
    out_node_id = out_node_id_; out_subnode_id = out_subnode_id_;
    import std.string : format;
    writeln("IN: <%d, %d>. OUT: <%d, %d>, ID: %d"
           .format(in_node_id_, in_subnode_id_, out_node_id_, out_subnode_id_, id_));
    writeln("RESULT: ", this.To_String);
  }

  this ( NodeConnection nc ) {
    in_node_id    = nc.in_node_id;       out_node_id = nc.out_node_id;
    in_subnode_id = nc.in_subnode_id; out_subnode_id = nc.out_subnode_id;
    id = nc.id;
  }

  auto To_String ( ) {
    import std.string : format;
    return "IN: <%d, %d>. OUT: <%d, %d>, ID: %d"
           .format(in_node_id, in_subnode_id, out_node_id, out_subnode_id, id);
  }
}


class UserValue {
  SubnodeDescription description;
  string value;
  this ( SubnodeDescription description_ ) {
    description = description_;
    switch ( description ) with ( SubnodeDescription ) {
      default: assert(0);
      case Colour: value = "[0.0, 0.0, 0.0]"; break;
      case Float3: value = "[0.0, 0.0, 0.0]"; break;
      case Float2: value = "[0.0, 0.0]";      break;
      case Float : value = "0.0";             break;
      case Int   : value = "0";               break;
      case String: value = "";                break;
    }
  }
}

public class Node {
private:
  ImVec2 origin, size;
  int id;
  string name;
  SubnodeData subnodes;
  UserValue user_value;

  void Calc_Subnode_Size(T)(T puts, ref float len) {
    import functional;
    len = 0.0f;
    foreach ( n; puts ) len = max(n.name.length, len);
    len *= Text_size;
  }

  void Calculate_Size ( ) {
    import functional : max, min;
    { // -- calculate width --
      // -- calculate input, output and title lengths in pixels --
      float title_length = name.length*Text_size,
            in_length    = RMax_Input_Width,
            out_length   = RMax_Output_Width;
      size.x = max(title_length, in_length + out_length);
    }
    { // -- calculate height --
      // -- calculate longest list from output/input --
      float list_len = max(subnodes.data.inputs.length,
                           subnodes.data.outputs.length);
      size.y = 25.0f + list_len*20.0f;
    }
    size.x = max(size.x + 45.0f,  85.0f);
    size.y = max(size.y +  5.0f,  25.0f);
    if ( user_value !is null ) {
      switch ( user_value.description ) with ( SubnodeDescription ) {
        default: assert(0);
        case Colour:
          size.x = 200.0f; break;
        case Float3: size.y = max(size.y, 100.0f); goto case;
        case Float2: size.y = max(size.y,  70.0f); goto case;
        case Int: case Float : size.y = max(size.y,  50.0f);
          size.x =  80.0f; break;
        case String: size.x = 100.0f; break;
      }
    }
  }

  auto RMax_Put_Size(T)(T put) {
    float len;
    Calc_Subnode_Size(put, len);
    return len;
  }
public:
  this ( NodeType node_type, int id_, ImVec2 origin_ ) {
    id = id_;
    subnodes = new SubnodeData(node_type);
    name = node_type.name;
    origin = origin_;
    if ( node_type.user_input != SubnodeDescription.Nil ) {
      user_value = new UserValue(node_type.user_input);
    }
    Calculate_Size();
  }

  this ( string node_type_name, int id_, ImVec2 origin_ ) {
    this(g_node_types[node_type_name], id_, origin_);
  }

  auto RID         ( ) inout      { return id;         }
  auto RName       ( ) inout      { return name;       }
  auto ROrigin     ( ) inout      { return origin;     }
  auto RSize       ( ) inout      { return size;       }
  auto RUser_Value ( )            { return user_value; }
  auto SUser_Value ( string val ) in {
    assert(user_value || val == "", "Setting val to null user value");
  } body {
    if ( user_value )
      user_value.value = val;
  }

  auto SOrigin ( ImVec2 origin_ ) {
    import std.math : fmax;
    origin = origin_;
  }

  auto RInput_Connection_Descs ( ) inout {
    return subnodes.type.subnode_descs.inputs;
  }
  auto ROutput_Connection_Descs ( ) inout {
    return subnodes.type.subnode_descs.outputs;
  }
  void SInput_Connection_Data  ( SubnodeDataContainer[] type ) {
    subnodes.data.inputs = type;
  }
  void SOutput_Connection_Data ( SubnodeDataContainer[] type ) {
    subnodes.data.outputs = type;
  }

  void Subnode_Connection ( int subnode_id, int connection_id, bool is_input ){
    (*(is_input ? &subnodes.data.inputs : &subnodes.data.outputs)
      )[subnode_id].data = connection_id;
  }

  // returns array tuple [bool is_input, SubnodeType]
  auto RSubnode_Connections ( ) {
    import functional, std.typecons : tuple;
    return subnodes.data.inputs .map!((ref n) => tuple(1, &n)).array ~
           subnodes.data.outputs.map!((ref n) => tuple(0, &n)).array;
  }
  auto RConnection_Descs ( ) {
    import functional, std.typecons : tuple;
    return RInput_Connection_Descs .map!(n => tuple(1, n)).array ~
           ROutput_Connection_Descs.map!(n => tuple(0, n)).array;
  }
  auto RSubnodes ( ) inout { return subnodes; }
  auto RMax_Output_Width ( ) { return RMax_Put_Size(subnodes.data.outputs); }
  auto RMax_Input_Width  ( ) { return RMax_Put_Size(subnodes.data.inputs ); }
}

public struct SubnodeDataContainer {
  string name;
  int data;
  this ( string name_, int connection_id_ ) {
    name = name_;
    data = connection_id_;
  }
}
class SubnodeData {
public:

  NodeType type;
  ConnectorBase!SubnodeDataContainer data;
  this ( NodeType type_ ) {
    import functional;
    type = type_;
    data.inputs = type_.subnode_descs.inputs.map!((SubnodeType n) {
      return SubnodeDataContainer(n.name, -1);
    }).array;
    data.outputs = type_.subnode_descs.outputs.map!((SubnodeType n) {
      return SubnodeDataContainer(n.name, -1);
    }).array;
    data.inputs.length  = type_.RInput_Length;
    data.outputs.length = type_.ROutput_Length;
  }
}

public enum SubnodeDescription {
  Colour, Float3, Float2, Float, Int, Varying, String,

  Texture, BRDF, Nil
};

struct SubnodeType {
  string             name;
  SubnodeDescription description;
  auto ROrigin ( ) {
    return ImVec2(name.length*13, 13);
  }

  bool RHovering ( ImVec2 offset ) {
    ImVec2 mouse;
    igGetMousePos(&mouse);
    auto p = Sub(mouse, Add(offset, ROrigin));
    writeln("P: ", p);
    return (p.x*p.x + p.y*p.y) < (Node_slot_radius*Node_slot_radius);
  }
}

public enum ParseNodeType {
  Function,
  Operation,
  Constant
}

class NodeType {
private:
  string name;
  ConnectorBase!SubnodeType subnode_descs;
  SubnodeDescription user_input;
  ParseNodeType parse_node_type;
public:
  this ( string name_, SubnodeType[] inputs, SubnodeType[] outputs,
         ParseNodeType parse_node_type_ ) {
    name = name_;
    subnode_descs = ConnectorBase!SubnodeType(inputs, outputs);
    user_input = SubnodeDescription.Nil;
    parse_node_type = parse_node_type_;
  }

  this ( string name_, SubnodeType[] inputs, SubnodeType[] outputs,
                       SubnodeDescription user_input_ ) {
    this(name_, inputs, outputs, ParseNodeType.Constant);
    user_input = user_input_;
  }

  auto RInput_Length    ( ) inout { return subnode_descs.inputs.length;  }
  auto ROutput_Length   ( ) inout { return subnode_descs.outputs.length; }
  auto RParse_Node_Type ( ) inout { return parse_node_type;              }
  string RName ( ) { return name; }
}
