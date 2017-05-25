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
    results ~= SubnodeDataContainer(load_subnode_data["name"].str, -1);
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
    auto load_nodes = load["nodes"];
    SNode_ID_Counter(load_nodes.array[0].str.to!int);
    foreach ( load_node; load_nodes.array[1 .. $] ) {
      auto node = New_Node(load_node["name"   ].str,
                           load_node["origin" ].str.To_ImVec2,
                           load_node["id"     ].str.to!int);
      node.SInput_Connection_Data (load_node["inputs" ].LPut_Data);
      node.SOutput_Connection_Data(load_node["outputs"].LPut_Data);
      node.SUser_Value(load_node["user_data"].str);
    }
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

string To_OpenCL_Float ( string float_val ) in {
  assert(!float_val.empty);
} body {
  import functional;
  // -- "3"|"3."|"3.0" -> "3.0f" (there are no f in imgui floats)
  if ( float_val.filter!(n => n == '.').empty ) float_val ~= ".";
  if ( float_val[$-1] == '.' ) float_val ~= "0"; // to include "3." & "3"
  return float_val ~ "f";
}

string Eval_Node_Str ( Node node ) {
  import functional;
  string nodename = node.RName;
  // -- special case nodename --
  switch ( nodename ) {
    default: break;
    case "Render_Map":
      return "res = opU(a, res, (float2)(%s, (float)(%s)));";
  }
  auto nodetype = nodename.RNodeType;
  final switch ( RNodeType(nodename).RParse_Node_Type ) {
    case ParseNodeType.Function:
      int amt = cast(int)nodetype.RInput_Length;
      auto arr = iota(0, amt).map!(n => "%s").array;
      if ( arr.length > 0 )
        arr[0..$-1] = arr[0..$-1].map!(n => n ~ ", ").array;
      return nodename ~ "(" ~ arr.joiner.to!string ~ ")";
    case ParseNodeType.Operation:
      return "(%s " ~ nodename ~ " %s)";
    case ParseNodeType.Constant:
      switch ( nodename ) {
        default: return nodename;
        case "Float": return node.RUser_Value.value.To_OpenCL_Float;
        case "Float2": case "Float3": case "Colour":
          auto vals = node.RUser_Value.value.to!(float[])
                          .map!(n => n.to!string.To_OpenCL_Float).array;
          if ( vals.length == 2 )
            return "(float2)(" ~ vals[0] ~ ", " ~ vals[1] ~ ")";
          return "(float3)(" ~ vals[0] ~ ", " ~ vals[1] ~ ", " ~ vals[2] ~ ")";
        case "String": case "Int": return node.RUser_Value.value;
        case "Origin": return "origin";
        case "Time":   return "time";
      }
  }
}

string Eval_Node ( Node node ) {
  string res = node.Eval_Node_Str;
  if ( node.RName.RNodeType.RParse_Node_Type == ParseNodeType.Constant ) {
    return res;
  }
  string[] args;
  foreach ( input_con; node.RSubnodes.data.inputs ) {
    auto out_node = input_con.data.RNodeConnection.out_node_id.RNode;
    auto ev = out_node.Eval_Node;
    args ~= ev;
  }
  if ( args.length == 1 ) res = res.format(args[0]);
  if ( args.length == 2 ) res = res.format(args[0], args[1]);
  if ( args.length == 3 ) res = res.format(args[0], args[1], args[2]);
  if ( args.length == 4 ) res = res.format(args[0], args[1], args[2], args[3]);
  return res;
}

string Parse_Graph ( ) {
  import functional;
  string parsed_result;
  // -- grab node map and gen str for each -- (TODO optimization, make O(1))
  foreach ( node; RNodes ) {
    if ( node.RName != "Render_Map" ) continue;
    parsed_result ~= "\n" ~ node.Eval_Node;
  }
  writeln("Result: ", parsed_result);
  import kernelinfo;
  Set_Map_Function(parsed_result);
  return "";
}
