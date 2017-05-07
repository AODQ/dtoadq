module gui.nodegraphrenderer;
import derelict.imgui.imgui;
import globals, gui.gui, gui.tnodes, gui.splines;
import DragState = gui.dragstate;
private:

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

void Update_Node ( ImDrawList* draw_list, ImVec2 offset, Node node,
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
  bool node_active = igIsItemActive() && DragState.Is_Drag_Default;
  auto node_bg_colour = node_hovered_in_scene == node.RID ?
                        ImColor(50,50,50) : ImColor(20,20,20);
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
  foreach ( it, ref connection; node.RConnection_Descs ) {
    auto subnode_text_size = gdCalcTextSize(connection.name);
    bool is_input = connection[0];
    auto subnode  = connection[1];
    auto pt = Add(Add(node_rect_min, connection.ROrigin),
                  is_input ?
                    ImVec2(connection.ROrigin.x - (text_size.x+1.0f), 0.0f) :
                    ImVec2(10.0f,                                    -6.0f));
    // -- render name and circle --
    igSetCursorScreenPos(pt);
    gdText(subnode.name);
    auto im_colour = ImColor(60, 60, 60);
    pt = Add(pt, ImVec2(is_input ? -10.0f : 20.0f, 0.0f));
    ImDrawList_AddCircleFilled(draw_list, pt, Node_slot_radius, im_colour);
    // -- check for input --
    import std.math : sqrt; // sdSphere { length(p) - radius }
    auto p = sqrt(pt.x*pt.x + pt.y*pt.y) - Node_slot_radius;
    if ( p <= 0.0f ) { DragState.Hover_Update(node.RID, it, is_input); }
    offset += 15.0f;
  }
  node_selected = node.RID;
  if (node_active && igIsMouseDragging(0) )
    node.SOrigin(Add(node.ROrigin, igGetIO().MouseDelta));
  // igEnd
  igPopId();
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

  // -- update nodes --
  foreach ( ref node; RNodes ) {
    Update_Node(draw_list, scrolling, node, node_selected);
    Render_Lines(draw_list, scrolling);
  }
  // -- update dragstate --
  DragState.Post_Update;

  // -- context menu --

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
