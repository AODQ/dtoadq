module gui.nodes;
import derelict.imgui.imgui;
import globals;
import gui.gui;
import std.variant : Variant;

/**
  Based heavily off of https://gist.github.com/emoon/b8ff4b4ce4f1b43e79f2
*/

immutable float Node_slot_radius = 5.0f;
immutable ImVec2 Node_window_padding = ImVec2(8.0f, 8.0f);
private DragNode drag_node;
private Node[] nodes;

auto Add ( inout ImVec2 x, inout ImVec2 y ) {
  return ImVec2(x.x + y.x, x.y + y.y);
}
auto Sub ( inout ImVec2 x, inout ImVec2 y ) {
  return ImVec2(x.x - y.x, x.y - y.y);
}

enum ConnectionType {
  Colour, Float3, Float2, Float, Int, Varying
}


struct ConnectionDesc {
  string name;
  ConnectionType type;
  this ( string name_ ) {
    name = name_;
    type = ConnectionType.Varying;
  }
}

// Can't use default constructor .. bug ? idk
private auto NConnectionDesc ( string name, ConnectionType type ) {
  ConnectionDesc connection_description;
  connection_description.name = name;
  connection_description.type = type;
  return connection_description;
}

struct NodeType {
  string name;
  ConnectionDesc[] input_connections;
  ConnectionDesc[] output_connections;
}

struct Varying {
  Variant data;
}

class Connection {
  ImVec2 origin;
  ConnectionDesc desc;
  Connection input_connection;
  Connection[] output_connections;

  Varying data;

  this ( ConnectionDesc desc_ ) {
    desc = desc_;
  }

  auto RInput ( ) { return input_connection; }

  auto ROutput ( size_t index ) in {
    assert(output_connections.length > index);
  } body {
    return output_connections[0];
  }
}

private alias _CD = ConnectionDesc, _CT = ConnectionType, _ND = NConnectionDesc;
private auto GNode ( string name, int input = 1 ) {
  import functional;
  return NodeType(name, iota(0, input).map!(n => _CD("In"~n.to!string)).array,
                        [_CD("Out")]);
}

private auto MNode ( string name, _CD[] arr ) {
  return NodeType(name, [_ND("Origin", _CT.Float3)] ~ arr,
                        [ _ND("Distance", _CT.Float) ]);
}

static immutable NodeType[] Node_types = [
  // -- basic generic operations --
  GNode("Multiply", 2), GNode("Add", 2),
  GNode("Divide", 2), GNode("Subtract", 2),
  NodeType("Random", [_CD("Low"), _CD("Hi")], [_CD("Output")]),
  GNode("Cos"), GNode("Sin"), GNode("Tan"),
  // -- map operations --
  MNode("sdSphere", [_ND("Radius", _CT.Float)]),
  MNode("sdBox",    [_ND("Size"  , _CT.Float)]),
  MNode("opRepeat", [_ND("Modulo", _CT.Float)]),
];


Connection[] Setup_Connections ( inout ConnectionDesc[] descs ) {
  import functional;
  return descs.dup.map!(n => new Connection(n)).array;
}

struct Node {
  ImVec2 origin, size;
  int id;
  string name;
  Connection[] input_connections;
  Connection[] output_connections;

  this ( int id_, string name_, ImVec2 origin_ ) {
    id = id_; name = name_; origin = origin_;
  }
}

Node Create_Node ( ImVec2 origin, inout NodeType node_type ) {
  static int id_counter = 0;
  Node node = Node(id_counter ++, node_type.name, origin);

  ImVec2 title_size = gdCalcTextSize(node.name);

  title_size.y *= 3;

  node.input_connections  = Setup_Connections(node_type.input_connections);
  node.output_connections = Setup_Connections(node_type.output_connections);

  auto input_text_size  = ImVec2(0.0f, 0.0f),
       output_text_size = input_text_size;
  foreach ( ref connection; node.input_connections ) {
    import std.math : fmax;
    ImVec2 text_size = gdCalcTextSize(node.name);
    input_text_size.x = fmax(text_size.x, input_text_size.x);

    connection.origin = ImVec2(0.0f, title_size.y + input_text_size.y +
                                     text_size.y/2.0f);

    input_text_size.y += text_size.y + 4.0f;
  }
  input_text_size.x +=  5.0f;


  float x_start = input_text_size.x;

  input_text_size.x = 0;
  foreach ( connection; node.input_connections ) {
    import std.math : fmax;
    auto text_size = gdCalcTextSize(connection.desc.name);
    input_text_size.x = fmax(text_size.x, input_text_size.x);
  }

  node.size.x = input_text_size.x + x_start;
  node.size.y = input_text_size.y + title_size.y;

  input_text_size.y = 0.0f;

  // set origin for output nodes if we know where to place

  foreach ( ref connection; node.output_connections ) {
    auto text_size = gdCalcTextSize(connection.desc.name);
    connection.origin = ImVec2(node.size.x, title_size.y + input_text_size.y
                                                         + text_size.y/2.0f);
    input_text_size.y += text_size.y + 4.0f;
  }

  return node;
}

auto New_Node ( ImVec2 origin, string name ) {
  nodes ~= Create_Node(origin, name);
}

auto Create_Node ( ImVec2 origin, string name ) {
  foreach ( ref node; Node_types ) {
    if ( node.name == name ) return Create_Node(origin, node);
  }
  assert(false);
}

auto Find_Node ( inout Connection connection ) {
  foreach ( ref node; nodes ) {
    foreach ( in_connection; node.input_connections ) {
      if ( in_connection == connection ) return &node;
    }
  }

  foreach ( ref node; nodes ) {
    foreach ( out_connection; node.output_connections ) {
      if ( out_connection == connection ) return &node;
    }
  }

  return null;
}

enum DragState {
  Default, Hover, Begin, Dragging, Connect
}

struct DragNode {
  ImVec2 origin;
  Connection connection;
  DragState state;
}


void Draw_Hermite ( ImDrawList* draw_list, ImVec2 p1, ImVec2 p2, int steps ) {
  ImVec2 t1 = ImVec2(80.0f, 0.0f), t2 = t1;

  foreach ( step; 0 .. steps ) {
    float t = float(step)/float(steps);
    float h1 =  2*t*t*t - 3*t*t + 1.0f,
          h2 = -2*t*t*t + 3*t*t,
          h3 =    t*t*t - 2*t*t + t,
          h4 =    t*t*t -   t*t;
    ImDrawList_PathLineTo(draw_list,
      ImVec2(h1*p1.x + h2*p2.x + h3*t1.x + h4*t2.x,
             h1*p1.y + h2*p2.y + h3*t1.x + h4*t2.y));
  }

  ImDrawList_PathStroke(draw_list, ImColor(150, 10, 175), false, 3.0f);
}


bool Is_Connector_Hovered ( inout ref Connection connection, ImVec2 offset ) {
  ImVec2 mouse_origin = gdRMousePos;
  ImVec2 connection_origin = Add(offset, connection.origin);

  float xd = mouse_origin.x - connection_origin.x,
        yd = mouse_origin.y - connection_origin.y;

  return ((xd*xd) + (yd*yd)) < (Node_slot_radius*Node_slot_radius);
}


Connection Get_Hover_Connection ( ImVec2 offset, out ImVec2 origin ) {
  foreach ( ref node; nodes ) {
    ImVec2 node_origin = Add(node.origin, offset);
    foreach ( ref connection; node.input_connections ) {
      if ( Is_Connector_Hovered(connection, node_origin) ) {
        origin = Add(node_origin, connection.origin);
        return connection;
      }
    }

    foreach ( ref connection; node.output_connections ) {
      if ( Is_Connector_Hovered(connection, node_origin) ) {
        origin = Add(node_origin, connection.origin);
        return connection;
      }
    }
  }

 return null;
}

void Update_Drag ( ImVec2 offset ) {
  ImVec2 origin;
  Connection connection = Get_Hover_Connection(offset, origin);
  final switch ( drag_node.state ) {
    case DragState.Default:
      if ( connection ) {
        drag_node.connection = connection;
        drag_node.origin     = origin;
        drag_node.state = DragState.Hover;
      }
    break;
    case DragState.Hover:
      if ( connection != drag_node.connection ) {
        drag_node.connection = null;
        drag_node.state = DragState.Default;
        return;
      }

      if ( igIsMouseClicked(0, 0) && drag_node.connection ) {
        drag_node.state = DragState.Dragging;
      }
    break;
    case DragState.Dragging:
      ImDrawList* draw_list = igGetWindowDrawList();
      ImDrawList_ChannelsSetCurrent(draw_list, 0); // set background
      Draw_Hermite(draw_list, drag_node.origin, gdRMousePos, 12);

      if ( !igIsMouseDown(0) ) {
        if ( connection == drag_node.connection ) {
          drag_node.connection = null;
          drag_node.state = DragState.Default;
          return;
        }

        drag_node.state = DragState.Default;
        if ( connection is null ) return;
        connection.input_connection = drag_node.connection;
        drag_node.connection = null;
        // Connect nodes (TODO: assure it's connected properly)
      }
    break;
    case DragState.Connect: case DragState.Begin: break;
  }

}

void Display_Node ( ImDrawList* draw_list, ImVec2 offset, ref Node node,
                    ref int node_selected ) {
  int node_hovered_in_scene = -1;
  bool open_context_menu = false;

  igPushIdInt(node.id);
    ImVec2 node_rect_min = Add(offset, node.origin);

    // Display node contents
    ImDrawList_ChannelsSetCurrent(draw_list, 1); // set foreground
    bool old_any_active = igIsAnyItemActive();

    // draw title in center
    ImVec2 text_size = gdCalcTextSize(node.name);
    ImVec2 origin = Add(node_rect_min, Node_window_padding);
    origin.x = node_rect_min.x + text_size.x / 2;
    igSetCursorScreenPos(origin);
    // igBeginGroup(); // Lock horizontal position
    gdText(node.name);

    bool node_widgets_active = (!old_any_active && igIsAnyItemActive());
    //node->size = ImGui::GetItemRectSize() + NODE_WINDOW_PADDING +
    //NODE_WINDOW_PADDING;
    auto node_rect_max = Add(node_rect_min, node.size);
    // Display node box
    ImDrawList_ChannelsSetCurrent(draw_list, 0); // set background
    igSetCursorScreenPos(node_rect_min);
    igInvisibleButton("node", node.size);
    if ( igIsItemHovered() ) {
      node_hovered_in_scene = node.id;
      open_context_menu |= igIsMouseClicked(1);
    }
    bool node_moving_active = (igIsItemActive() && !drag_node.connection);
    auto node_bg_colour = node_hovered_in_scene == node.id ?
                          ImColor(50,50,50) :
                          ImColor(20,20,20);
    ImDrawList_AddRectFilled(draw_list, node_rect_min, node_rect_max,
                             node_bg_colour, 4.0f, 0);
    auto title_area = node_rect_max;
    title_area.y = node_rect_min.y + 3.0f;
  // Draw text bg area
    ImDrawList_AddRectFilled(draw_list, Add(node_rect_min, ImVec2(1,1)),
                                 title_area, ImColor(100,0,0), 4.0f, 0);
    ImDrawList_AddRect(draw_list, node_rect_min, node_rect_max,
                                 ImColor(50, 50, 50), 4.0f, 0, 1.0f);
    ImVec2 off = node_rect_min;
    offset.y += 40.0f;
    offset = Add(offset, node_rect_min);
    foreach ( ref connection; node.input_connections ) {
      auto res = Add(Add(node_rect_min, connection.origin),
                     ImVec2(10.0f, -6.0f));
      igSetCursorScreenPos(res);
      gdText(connection.desc.name);
      auto im_colour = ImColor(60, 60, 60);
      if (Is_Connector_Hovered(connection, node_rect_min))
        im_colour = ImColor(80, 80, 80);
        ImDrawList_AddCircleFilled(draw_list,
          Add(node_rect_min, connection.origin), Node_slot_radius, im_colour);
      offset.y += text_size.y + 2.0f;
    }

    offset = node_rect_min;
    offset.y += 40.0f;
    foreach ( ref connection; node.output_connections ) {
      text_size = gdCalcTextSize(connection.desc.name);
      igSetCursorScreenPos(Add(offset,
        ImVec2(connection.origin.x - (text_size.x + 1.0f), 0)));
      gdText(connection.desc.name);
      auto im_colour = ImColor(60, 60, 60);
      if (Is_Connector_Hovered(connection, node_rect_min))
        im_colour = ImColor(80, 80, 80);
      ImDrawList_AddCircleFilled(draw_list,Add(node_rect_min,connection.origin),
                                Node_slot_radius, im_colour);
      offset.y += text_size.y + 2.0f;
    }
    if (node_widgets_active || node_moving_active)
      node_selected = node.id;
    if (node_moving_active && igIsMouseDragging(0) )
      node.origin = Add(node.origin, igGetIO().MouseDelta);
    // igEndGroup();
  igPopId();
}

void Render_Lines ( ImDrawList* draw_list, ImVec2 offset ) {
  foreach ( ref node; nodes ) {
    foreach ( ref connection; node.input_connections ) {
      if ( connection.RInput is null ) continue;
      auto target_node = Find_Node(connection.RInput);
      if ( !target_node ) continue;

      Draw_Hermite(draw_list,
        Add(Add(offset, target_node.origin), connection.RInput.origin),
        Add(Add(offset, node.origin),        connection.origin), 12);
    }
  }
}



void Show_Example_Node_Graph ( bool* opened ) {
  igSetNextWindowSize(ImVec2(700, 600), ImGuiSetCond_FirstUseEver);
  if ( !igBegin("Example: Custom node graph", opened) ) {
    igEnd();
    return;
  }

  bool open_context_menu = false;
  int node_hovered_in_list  = -1,
      node_hovered_in_scene = -1;
  int node_selected = -1;
  static ImVec2 scrolling = ImVec2(0.0f, 0.0f);

  igSameLine();
  igBeginGroup();

  igPushStyleVarVec(ImGuiStyleVar_FramePadding, ImVec2(1, 1));
  igPushStyleVarVec(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
  // igPushStyleColor(ImGuiCol_ChildWindowBg, ImVec4(255, 20, 20, 200));
  igBeginChild("scrolling_region", ImVec2(0, 0), true,
                ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoMove);
  igPushItemWidth(120.0f);
  auto draw_list = igGetWindowDrawList();
  ImDrawList_ChannelsSplit(draw_list, 2);

  foreach ( ref node; nodes ) {
    Display_Node(draw_list, scrolling, node, node_selected);
    Update_Drag(scrolling);
    Render_Lines(draw_list, scrolling);
  }

  static ImVec2 saved_mouse_pos;
  // open context menu
  if ( igIsAnyItemHovered && igIsMouseHoveringWindow && igIsMouseClicked(1)){
    node_selected = node_hovered_in_list = node_hovered_in_scene = -1;
    open_context_menu = true;
    saved_mouse_pos = gdRMousePos;
  }
  if ( open_context_menu ) {
    igOpenPopup("context_menu");
    if ( node_hovered_in_list  != -1 ) node_selected = node_hovered_in_list;
    if ( node_hovered_in_scene != -1 ) node_selected = node_hovered_in_scene;
  }
  // draw context menu
  igPushStyleVarVec(ImGuiStyleVar_WindowPadding, ImVec2(8, 8));
    if ( igBeginPopup("context_menu") ) {
      if ( gdMenuItem("Load Graph") ) {
        // TODO
      }
      if ( gdMenuItem("Save Graph") ) {
        writeln("Saved graph: ");
        string results;
        foreach ( node; nodes ) {
          results ~= "NODE ---\n";
          results ~= "  Origin: <" ~ node.origin.x.to!string ~ ", "
                                 ~ node.origin.y.to!string ~ ">,\n";
          results ~= "  Name: " ~ node.name ~ "\n";
          results ~= "  Input: TODO\n";
        }
        writeln(results);
      }

      if ( igBeginMenu("New Node") ) {
        if ( igBeginMenu("Generic Math") ) {
          if ( gdMenuItem("Multiply") ) {
            New_Node(saved_mouse_pos, "Multiply");
          }
          if ( gdMenuItem("Cos") ) {
            New_Node(saved_mouse_pos, "Cos");
          }
          if ( gdMenuItem("Sin") ) {
            New_Node(saved_mouse_pos, "Sin");
          }
          if ( gdMenuItem("Tan") ) {
            New_Node(saved_mouse_pos, "Tan");
          }
          if ( gdMenuItem("Subtract") ) {
            New_Node(saved_mouse_pos, "Subtract");
          }
          if ( gdMenuItem("Divide") ) {
            New_Node(saved_mouse_pos, "Divide");
          }
          if ( gdMenuItem("Add") ) {
            New_Node(saved_mouse_pos, "Add");
          }
          igEndMenu();
        }
        if ( igBeginMenu("Map Operations") ) {
          if ( gdMenuItem("sdSphere") ) {
            New_Node(saved_mouse_pos, "sdSphere");
          }
          if ( gdMenuItem("sdBox") ) {
            New_Node(saved_mouse_pos, "sdBox");
          }
          if ( gdMenuItem("opRepeat") ) {
            New_Node(saved_mouse_pos, "opRepeat");
          }
          igEndMenu();
        }
        igEndMenu();
      }

      igEndPopup();
    }
  igPopStyleVar();

  ImDrawList_ChannelsMerge(draw_list);

  // Scrolling
  if ( igIsWindowHovered() && !igIsAnyItemActive()&&igIsMouseDragging(2, 0.0f))
    scrolling = Sub(scrolling, igGetIO().MouseDelta);

  igPopItemWidth();
  igEndChild();
  // igPopStyleColor();
  igPopStyleVar();
  igPopStyleVar();
  igEndGroup();
  igEnd();
}
