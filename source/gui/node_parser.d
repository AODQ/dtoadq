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
  JSONValue save_subnode_puts = JSONValue([puts.length.to!string]);
  foreach ( ref subnode_data; puts ) {
    JSONValue save_subnode_data = [
      "name" : subnode_data.name.to!string, // not 'necessary' but easy for now
      "data" : subnode_data.data.to!string,
    ];
    save_subnode_puts.array ~= save_subnode_data;
  }
  return save_subnode_puts;
}

auto LPut_Data(T)( T load_data ) {
  import std.json;
  auto len = load_data.array[0].str.to!int;
  SubnodeDataContainer[] results;
  foreach ( load_subnode_data; load_data.array[1 .. $] ) {
    results ~= SubnodeDataContainer(load_subnode_data["name"].str,
                                    load_subnode_data["data"].str.to!int);
  }
  return results;
}

void Save_Graph ( ) {
  import std.json;
  // I convert everything to strings to avoid headaches from converting from
  // floats to ints etc, I just remember it always being a hassle otherwise
  JSONValue save = ["Version" : "0.0"];
  // Have to save nodes & node connections
  { // -- save nodes --
    save.object["nodes"] = JSONValue([RNode_ID_Counter.to!string]);
    foreach ( node; RNodes ) {
      auto subnodes = node.RSubnodes;
      JSONValue save_node = [
        "id"      : node.RID.to!string,
        "origin"  : node.ROrigin.To_String,
        "name"    : node.RName,
      ];
      save_node.object["inputs" ] = RPut_Data(subnodes.data.inputs);
      save_node.object["outputs"] = RPut_Data(subnodes.data.outputs);
      // -- user data --
      auto user_val = node.RUser_Value;
      writeln("User value: ", user_val);
      save_node.object["user_data"] = "";
      if ( user_val ) {
        save_node["user_data"] = user_val.value;
      }
      save["nodes"].array ~= save_node;
    }
  }

  { // -- save connections --
    save.object["connections"]=JSONValue([RNode_Connections_Counter.to!string]);
    foreach ( con; RNode_Connections ) {
      save["connections"].array ~= JSONValue([
        "id"             : con.id             .to!string,
        "in_node_id"     : con.in_node_id     .to!string,
        "out_node_id"    : con.out_node_id    .to!string,
        "in_subnode_id"  : con.in_subnode_id  .to!string,
        "out_subnode_id" : con.out_subnode_id .to!string,
      ]);
    }
  }

  static import std.file;
  std.file.write("temp-node-graph.json", save.toString);
}

void Clear_Graph ( ) {
  Clear_Nodes();
  Clear_Node_Connections();
}

void Load_Graph ( ) {
  Clear_Graph();
  import std.json;
  static import std.file;
  JSONValue load = parseJSON(std.file.read("temp-node-graph.json").to!string);

  { // -- load nodes --
    writeln("LOAD NODE");
    auto load_nodes = load["nodes"];
    writeln("SETING ID: ", load_nodes);
    SNode_ID_Counter(load_nodes.array[0].str.to!int);
    writeln("NOIDE ID COUNTER SET");
    foreach ( load_node; load_nodes.array[1 .. $] ) {
      writeln("NEW NODE");
      auto node = New_Node(load_node["name"   ].str,
                           load_node["origin" ].str.To_ImVec2,
                           load_node["id"     ].str.to!int);
      writeln("RESULTS: ", node);
      node.SInput_Connection_Data (load_node["inputs" ].LPut_Data);
      writeln("INPUT CONNECTION");
      node.SOutput_Connection_Data(load_node["outputs"].LPut_Data);
      writeln("OUTPUT CONNECTION");
      writeln("DATA: ", load_node["user_data"]);
      node.SUser_Value(load_node["user_data"].str);
      writeln("USER VALUE");
    }
    load["nodes"].writeln;
  }

  { // -- load connections --
    auto load_cons = load["connections"];
    SNode_Connections_Counter(load_cons.array[0].str.to!int);
    foreach ( load_con; load_cons.array[1..$] ) {
      auto id = load_con["id"].str.to!int;
      Add_Connection(new NodeConnection(
        id,
        load_con["in_node_id"     ].str.to!int,
        load_con["out_node_id"    ].str.to!int,
        load_con["in_subnode_id"  ].str.to!int,
        load_con["out_subnode_id" ].str.to!int), id);
    }
  }
}


struct GraphParseInfo {
  bool compiled;
  string result; // result if compiled, otherwise error string
}

string Parse_Graph ( ) {
  return "";
}
