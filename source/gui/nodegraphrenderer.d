module gui.nodegraphrenderer;
import derelict.imgui.imgui;
import globals, gui.gui, gui.tnodes, gui.splines;
import DragState = gui.dragstate;
private:

void Render_Lines ( ImDrawList* draw_list, ImVec2 offset ) {
  foreach ( ref con; RNode_Connections ) {
    auto in_node  = RNode(con.in_node_id),
         out_node = RNode(con.out_node_id);
    auto in_offset  = ImVec2(-8.0f, 25.0f + con.in_subnode_id*20.0f + 6.0f),
         out_offset = ImVec2(out_node.RSize.x + 8.0f,
                             25.0f + con.out_subnode_id*20.0f + 6.0f);
    Draw_Hermite(draw_list,
      Add(Add(offset, out_node.ROrigin), out_offset),
      Add(Add(offset, in_node .ROrigin),  in_offset)
    );
  }
}

void Update_Node ( ImDrawList* draw_list, ImVec2 screen_offset, Node node,
                   ref int node_selected ) {
  // -- mouse/ui info --
  bool old_any_active = igIsAnyItemActive();
  int node_hovered = -1;
  // -- calculate node bounds --
  auto node_bound_min = Add(screen_offset,  node.ROrigin),
       node_bound_max = Add(node_bound_min, node.RSize);

  igPushIdInt(cast(int)node.RID);

  auto user_value = node.RUser_Value;
  if ( user_value ) { // -- draw user input --
    // Have to draw this before noede bounds as node bounds uses button that
    // would override this input
    ImDrawList_ChannelsSetCurrent(draw_list, 1);
    auto pt = Add(node_bound_min, ImVec2(5.0f, 25.0f));
    igSetCursorScreenPos(pt);
    switch ( user_value.description ) with ( SubnodeDescription ) {
      default: assert(0);
      case Colour:
        float[3] arr = user_value.value.to!(float[]);
        igColorEdit3("C", arr);
        user_value.value = arr.to!string;
      break;
      case Float3:
        float[3] arr = user_value.value.to!(float[]);
        gdInputFloat("X", arr[0]);
        igSetCursorScreenPos(Add(pt, ImVec2(0.0f, 25.0f)));
        gdInputFloat("Y", arr[1]);
        igSetCursorScreenPos(Add(pt, ImVec2(0.0f, 50.0f)));
        gdInputFloat("Z", arr[2]);
        user_value.value = arr.to!string;
      break;
      case Float2:
        float[2] arr = user_value.value.to!(float[]);
        gdInputFloat("X", arr[0]);
        igSetCursorScreenPos(Add(pt, ImVec2(0.0f, 25.0f)));
        gdInputFloat("Y", arr[1]);
        user_value.value = arr.to!string;
      break;
      case Float :
        float val = user_value.value.to!float;
        gdInputFloat("F", val);
        user_value.value = val.to!string;
      break;
      case Int   :
        int val = user_value.value.to!int;
        gdInputInt("I", val);
        user_value.value = val.to!string;
      break;
      case String:
        // gdInputText("Text", user_value.value);
      break;
    }
  }

  { // -- draw node bounds --
    ImDrawList_ChannelsSetCurrent(draw_list, 0);
    igSetCursorScreenPos(node_bound_min);
    auto chr = "node " ~ node.RName;
    igInvisibleButton(chr.toStringz, node.RSize);
    if ( igIsItemHovered() && DragState.RDrag_State == DragState.State.Default){
      node_hovered      = cast(int)node.RID;
      node_selected     = igIsMouseClicked(1) ? node_hovered : node_selected;
    }
    auto bg_col = node_hovered >= 0 ? 0.2f : 0.1f;
    auto node_bg_colour = ImColor(bg_col, bg_col, bg_col);
    ImDrawList_AddRectFilled(draw_list, node_bound_min, node_bound_max,
                             node_bg_colour, 4.0f, 0);
  }

  ImDrawList_ChannelsSetCurrent(draw_list, 1);

  { // -- draw centered title --
    auto text_size = node.RName.length*Text_size;
    auto origin    = Add(node_bound_min, ImVec2(2.0f, 2.0f));
    origin.x += node.RSize.x/2.0f; // center
    origin.x -= text_size   /2.0f; // left adjust
    igSetCursorScreenPos(origin);
    gdText(node.RName);
  }

  { // -- draw input/outputs --
    auto coffset = ImVec2(3.0f, 25.0f);
    bool is_input = true;
    int last_input = 0;
    foreach ( it, ref connection; node.RConnection_Descs ) {
      if ( connection[0] != is_input ) {
        coffset.y =  25.0f;
        coffset.x = node.RSize.x - node.RMax_Output_Width - 5.0f;
        last_input = cast(int)it;
      }
      is_input     = cast(bool)connection[0];
      auto subnode = connection[1];
      auto subnode_text_size = subnode.name.length*Text_size;
      { // -- render name --
        auto pt = Add(node_bound_min, coffset);
        igSetCursorScreenPos(pt);
        gdText(subnode.name);
      }
      { // -- render circle --
        auto pt = Add(node_bound_min, ImVec2(0.0f, coffset.y + 6.0f));
        pt.x += is_input ? -8.0f : node.RSize.x + 8.0f;
        import std.math : sqrt; // I can use SDF of a sphere here
        auto m = Sub(pt, gdRMousePos);
        static Radius = 5.0f;
        bool hovering = (sqrt(m.x*m.x + m.y*m.y) - Radius) <= 0.0f;
        auto col = hovering ? 0.8f : 0.2f;
        auto colour = ImColor(col, col, col);
        ImDrawList_AddCircleFilled(draw_list, pt, Radius, colour);
        // -- update circle --
        if ( hovering )
          DragState.Hover_Update(node.RID, cast(int)it - last_input, is_input);
      }
      coffset.y += 20.0f;
    }
  }

  // -- check node movement --
  if ( node_hovered >= 0 && igIsMouseDragging(0) &&
       DragState.RDrag_State == DragState.DragState.Default )
    node.SOrigin(Add(node.ROrigin, igGetIO().MouseDelta));

  igPopId();
}

public void Update_Node_Graph ( ) {
  static bool opened = true;
  igSetNextWindowSize(ImVec2(400, 400), ImGuiSetCond_FirstUseEver);
  if ( !igBegin("SDF Node Graph", &opened) ) {
    igEnd();
    return;
  }

  bool open_context_menu = false;
  int node_hovered_in_list  = -1,
      node_hovered_in_scene = -1;
  static int node_selected = -1;
  static ImVec2 scrolling = ImVec2(0.0f, 0.0f);

  igSameLine();
  igBeginGroup();

  igPushStyleVarVec(ImGuiStyleVar_FramePadding,  ImVec2(1, 1));
  igPushStyleVarVec(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
  igBeginChild("scrolling_region", ImVec2(0, 0), true,
                ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoMove);
  igPushItemWidth(60.0f); // set input field
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
  if ( igIsMouseHoveringWindow && igIsMouseClicked(1)){
    open_context_menu = true;
    saved_mouse_pos = gdRMousePos;
  }
  if ( open_context_menu ) {
    if ( node_selected != -1 ) igOpenPopup("node_menu"   );
    else                       igOpenPopup("context_menu");
  }
  // draw context menu
  igPushStyleVarVec(ImGuiStyleVar_WindowPadding, ImVec2(8, 8));
    if ( node_selected != -1 ) {
      if ( igBeginPopup("node_menu") ) {
        if ( igBeginMenu("Erase") ) {
          if ( gdMenuItem("Node") ) {
            Remove_Node(node_selected);
          }
          if ( gdMenuItem("Connections") ) {
            Clear_Connections(node_selected);
          }
          igEndMenu();
        }
        igEndPopup();
      } else {
        node_selected = -1;
      }
    } else {
      if ( igBeginPopup("context_menu") ) {
        if ( gdMenuItem("Save Graph") ) {
          import gui.node_parser;
          Save_Graph();
        }
        if ( gdMenuItem("Load Graph") ) {
          import gui.node_parser;
          Load_Graph();
        }
        if ( gdMenuItem("Parse Graph") ) {
          import gui.node_parser;
          Parse_Graph();
        }

        void Menu ( string label, string[] names ) {
          if ( !igBeginMenu(label.toStringz) ) return;
          foreach ( name; names )
            if ( gdMenuItem(name) )
              New_Node(name, Sub(saved_mouse_pos, scrolling));
          igEndMenu();
        }

        if ( igBeginMenu("New Node") ) {
          import gui.opencl_funcs;
          Apply_Function(&Menu);
          igEndMenu();
        }

        igEndPopup();
      }
    }
  igPopStyleVar();

  ImDrawList_ChannelsMerge(draw_list);

  // Scrolling
  if ( igIsWindowHovered() && !igIsAnyItemActive()&&igIsMouseDragging(2, 0.0f))
    scrolling = Add(scrolling, igGetIO().MouseDelta);

  igPopItemWidth();
  igEndChild();
  igPopStyleVar();
  igPopStyleVar();
  igEndGroup();
  igEnd();
}
