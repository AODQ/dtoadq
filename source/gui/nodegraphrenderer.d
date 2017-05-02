module gui.nodegraphrenderer;
import derelict.imgui.imgui;
import globals, gui.gui;
import gui.tnodes;
private:

enum DragState {
  Default, Hover, Begin, Dragging, Connect
}
private ImVec2 g_drag_origin;
private NodeConnection g_drag_connection;
private DragState g_drag_state;

void Draw_Hermite ( ImDrawList* draw_list, ImVec2 p1, ImVec2 p2 ) {
  ImVec2 v1 = ImVec2(80.0f, 0.0f), v2 = v1;
  foreach ( step; 0 .. 9 ) {
    float t = cast(float)step/9.0f, t2 = t*t, t3 = t*t*t;
    float h1 =  2*t3 - 3*t2 + 1.0f,
          h2 = -2*t3 + 3*t2,
          h3 =    t3 - 2*t2 + t,
          h4 =    t3 -   t2;
    ImDrawList_PathLineTo(draw_list,
      ImVec2(h1*p1.x + h2*p2.x + h3*v1.x + h4*v2.x,
             h1*p1.y + h2*p2.y + h3*v1.x + h4*v2.y));
  }
  ImDrawList_PathStroke(draw_list, ImColor(150, 10, 175), false, 3.0f);
}

void Render_Lines ( ImDrawList* draw_list, ImVec2 offset ) {
  foreach ( ref node; RNodes ) {
    foreach ( it, ref subnode; node.RSubnode_Data.data.inputs ) {
      // if ( subnode is null ) continue;
      if ( subnode.type != SubnodeDataType.Connection ) continue;
      auto connection = RNodeConnection(subnode.data.connection_id);
      auto out_node = RNode(connection.out_node_id),
           out_subnode_id = connection.out_subnode_id;
      Draw_Hermite(draw_list,
        Add(Add(offset, out_node.ROrigin), ImVec2(           2, 15+it*10)),
        Add(Add(offset,     node.ROrigin), ImVec2(node.RSize.x, 15+it*10)));
    }
  }
}

void Render_Node ( ImDrawList* draw_list, ImVec2 offset, Node node,
                   ref int node_selected ) {
  int node_hovered_in_scene = -1;
  bool open_context_menu = false;

  igPushIdInt(cast(int)node.RID);
    ImVec2 node_rect_min = Add(offset, node.ROrigin);

    // Display node contents
    ImDrawList_ChannelsSetCurrent(draw_list, 1); // set foreground
    bool old_any_active = igIsAnyItemActive();

    // draw title in center
    ImVec2 text_size = gdCalcTextSize(node.RName);
    ImVec2 origin = Add(node_rect_min, Node_window_padding);
    origin.x = node_rect_min.x + text_size.x / 2;
    igSetCursorScreenPos(origin);
    gdText(node.RName);

    bool node_widgets_active = (!old_any_active && igIsAnyItemActive());
    auto node_rect_max = Add(node_rect_min, node.RSize);
    // Display node box
    ImDrawList_ChannelsSetCurrent(draw_list, 0); // set background
    igSetCursorScreenPos(node_rect_min);
    igInvisibleButton("node", node.RSize);
    if ( igIsItemHovered() ) {
      node_hovered_in_scene = cast(int)node.RID;
      open_context_menu |= igIsMouseClicked(1);
    }
    bool node_moving_active = (igIsItemActive() && !g_drag_connection);
    auto node_bg_colour = node_hovered_in_scene == node.RID ?
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
    foreach ( ref connection; node.RInput_Connection_Descs ) {
      auto res = Add(Add(node_rect_min, connection.ROrigin),
                     ImVec2(10.0f, -6.0f));
      igSetCursorScreenPos(res);
      gdText(connection.name);
      auto im_colour = ImColor(60, 60, 60);
      // if (Is_Connector_Hovered(connection, node_rect_min))
      //   im_colour = ImColor(80, 80, 80);
      ImDrawList_AddCircleFilled(draw_list,
          Sub(node_rect_min, connection.ROrigin), Node_slot_radius, im_colour);
      offset.y += text_size.y + 2.0f;
    }

    offset = node_rect_min;
    offset.y += 40.0f;
    foreach ( it, ref connection; node.ROutput_Connection_Descs ) {
      text_size = gdCalcTextSize(connection.name);
      igSetCursorScreenPos(Add(offset,
        ImVec2(connection.ROrigin.x - (text_size.x + 1.0f), 0)));
      gdText(connection.name);
      auto im_colour = ImColor(60, 60, 60);
      // if (Is_Connector_Hovered(connection, node_rect_min))
      //   im_colour = ImColor(80, 80, 80);
      ImDrawList_AddCircleFilled(draw_list,
           Add(node_rect_max, ImVec2(connection.ROrigin.x, it*10.0f)),
           Node_slot_radius, im_colour);
      offset.y += text_size.y + 2.0f;
    }
    if (node_widgets_active || node_moving_active)
      node_selected = node.RID;
    if (node_moving_active && igIsMouseDragging(0) )
      node.SOrigin(Add(node.ROrigin, igGetIO().MouseDelta));
    // igE
  igPopId();
}

void Update_Input_Drag ( ImVec2 offset ) {
  ImVec2 origin;
  auto node_connection = RHovering_Node(offset, origin);
  final switch ( g_drag_state ) {
    case DragState.Default:
      if ( g_drag_connection ) {
        g_drag_connection = node_connection;
        g_drag_origin     = origin;
        g_drag_state = DragState.Hover;
      }
    break;
    case DragState.Hover:
      if ( node_connection.in_node_id    != g_drag_connection.in_node_id ||
           node_connection.in_subnode_id != g_drag_connection.in_subnode_id ) {
        g_drag_connection = null;
        g_drag_state = DragState.Default;
        return;
      }

      if ( igIsMouseClicked(0, 0) && g_drag_connection )
        g_drag_state = DragState.Dragging;
    break;
    case DragState.Dragging:
      ImDrawList* draw_list = igGetWindowDrawList();
      ImDrawList_ChannelsSetCurrent(draw_list, 0); // set background
      Draw_Hermite(draw_list, g_drag_origin, gdRMousePos);

      if (!igIsMouseDown(0) ) break;
      g_drag_state = DragState.Default;
      if ( node_connection.in_node_id    == g_drag_connection.in_node_id &&
           node_connection.in_subnode_id == g_drag_connection.in_subnode_id ) {
        g_drag_connection = null;
        return;
      }

      // g_drag_connection
      // connection.input_connection = drag_node.connection;
      // drag_node.connection = null;
    break;
    case DragState.Connect: case DragState.Begin: break;
  }

}


public void Update_Node_Graph ( ) {
  static bool opened = true;
  igSetNextWindowSize(ImVec2(400, 400), ImGuiSetCond_FirstUseEver);
  if ( !igBegin("Model Node Graph", &opened) ) {
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
  igBeginChild("scrolling_region", ImVec2(0, 0), true,
                ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoMove);
  igPushItemWidth(120.0f);
  auto draw_list = igGetWindowDrawList();
  ImDrawList_ChannelsSplit(draw_list, 2);

  foreach ( ref node; RNodes ) {
    Render_Node(draw_list, scrolling, node, node_selected);
    Update_Input_Drag(scrolling);
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
        foreach ( node; RNodes ) {
          // results ~= "NODE ---\n";
          // results ~= "  Origin: <" ~ node.origin.x.to!string ~ ", "
          //                        ~ node.origin.y.to!string ~ ">,\n";
          // results ~= "  Name: " ~ node.name ~ "\n";
          // results ~= "  Input: TODO\n";
        }
        writeln(results);
      }

      if ( igBeginMenu("New Node") ) {
        if ( igBeginMenu("Generic Math") ) {
          if ( gdMenuItem("Multiply") ) {
            New_Node("Multiply", saved_mouse_pos);
          }
          if ( gdMenuItem("Cos") ) {
            New_Node("Cos", saved_mouse_pos);
          }
          if ( gdMenuItem("Sin") ) {
            New_Node("Sin", saved_mouse_pos);
          }
          if ( gdMenuItem("Tan") ) {
            New_Node("Tan", saved_mouse_pos);
          }
          if ( gdMenuItem("Subtract") ) {
            New_Node("Subtract", saved_mouse_pos);
          }
          if ( gdMenuItem("Divide") ) {
            New_Node("Divide", saved_mouse_pos);
          }
          if ( gdMenuItem("Add") ) {
            New_Node("Add", saved_mouse_pos);
          }
          igEndMenu();
        }
        if ( igBeginMenu("Map Operations") ) {
          if ( gdMenuItem("sdSphere") ) {
            New_Node("sdSphere", saved_mouse_pos);
          }
          if ( gdMenuItem("sdBox") ) {
            New_Node("sdBox", saved_mouse_pos);
          }
          if ( gdMenuItem("opRepeat") ) {
            New_Node("opRepeat", saved_mouse_pos);
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
