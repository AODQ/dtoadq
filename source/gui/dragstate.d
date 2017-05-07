module gui.dragstate;
import derelict.imgui.imgui, gui.tnodes, gui.splines, gui.gui;

enum DragState {
  Default, Hover, Dragging
}
private ImVec2 g_origin;
private NodeConnection g_connection;
private DragState g_state;
private bool g_is_input;
private bool last_update_hovered;

bool Valid_Subnode_Connection ( int node_id, int subnode_id, bool is_input ) {
  return is_input != g_is_input && g_connection.RNode_ID(is_input) != node_id &&
         g_connection.RSubnode_ID(is_input) != subnode_id;
}

void Reset_Connection ( ) {
  g_state      = DragState.Default;
  g_connection = null;
}

void Hover_Update ( int node_id, int subnode_id, bool is_input ) {
  final switch ( g_state ) with ( DragState ) {
    case Default:
      if ( !igIsMouseClicked(0) ) g_state = DragState.Hover;
    break;
    case Hover:
      if ( igIsMouseClicked(0) ) {
        g_state      = DragState.Dragging;
        g_connection = new NodeConnection(node_id, subnode_id, is_input);
        g_is_input   = is_input;
        g_origin     = gdRMousePos;
        last_update_hovered = true;
      }
    break;
    case Dragging:
      if ( !igIsMouseClicked(0) ) {
        if ( !Valid_Subnode_Connection(node_id, subnode_id, is_input) ) {
          Add_Connection(g_connection);
        }
        Reset_Connection;
      }
    break;
  }
}

void Post_Update ( ) {
  final switch ( g_state ) with ( DragState ) {
    case Default: break;
    case Hover:
      if ( !last_update_hovered ) g_state = DragState.Default;
    break;
    case Dragging:
      if ( !igIsMouseClicked(0) ) {
        Reset_Connection;
      }
      ImDrawList* draw_list = igGetWindowDrawList();
      ImDrawList_ChannelsSetCurrent(draw_list, 0); // set background
      Draw_Hermite(draw_list, g_origin, gdRMousePos);
    break;
  }
  last_update_hovered = false;
}


auto RDrag_State ( ) { return g_state; }
auto Is_Drag_Default ( ) { return g_state == DragState.Default; }
