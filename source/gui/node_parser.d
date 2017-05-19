module gui.node_parser;
import globals, gui.tnodes, derelict.imgui.imgui : ImVec2;

string To_String ( ImVec2 vec ) {
  return [vec.x, vec.y].to!string;
}

ImVec2 To_ImVec2 ( string val ) {
  auto arr = val.to!(float[]);
  return ImVec2(arr[0], arr[1]);
}

auto RPut_Data(T)(T puts) {
  import std.json : JSONValue;
  JSONValue[] save_subnode_puts;
  foreach ( immutable ref subnode_data; puts ) {
    JSONValue save_subnode_data;
    save_subnode_data.object["type"] = subnode_data.type.to!string;
    save_subnode_data.object["name"] = subnode_data.name.to!string;
    save_subnode_data.object["data"] = subnode_data.data.to!string;
    save_subnode_puts ~= save_subnode_data;
  }
  return save_subnode_puts;
}

void Save_Graph ( ) {
  import std.json;
  // I convert everything to strings to avoid headaches from converting from
  // floats to ints etc, I just remember it always being a hassle otherwise
  JSONValue save;
  // Have to save nodes & node connections
  { // save nodes
    save.object["node_id_counter"] = RNode_ID_Counter;
    foreach ( node; RNodes ) {
      JSONValue save_node;
      auto subnodes = node.RSubnodes;
      save_node.object["id"]     = node.RID.to!string;
      save_node.object["origin"] = node.ROrigin.To_String;
      save_node.object["name"] = node.RName;
      save_node.object["type"] = subnodes.type.to!string;
      save_node.object["inputs"] = RPut_Data(subnodes.data.inputs);
      save_node.object["outputs"] = RPut_Data(subnodes.data.outputs);
      { // -- user data --
        JSONValue save_user_val;
        auto user_val = node.RUser_Value;
        save_user_val.object["description"] = user_val.description.to!string;
        save_user_val.object["value"] = user_val.value; // already string
      }
    }
  }
}

void Load_Graph ( ) {
}


struct GraphParseInfo {
  bool compiled;
  string result; // result if compiled, otherwise error string
}

string Parse_Graph ( ) {
  return "";
}
