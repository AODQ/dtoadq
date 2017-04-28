module gui.tnodes;
import derelict.imgui.imgui;
import globals, gui.gui, std.variant : Variant;

class NodeTypeParameter {
  string name;
  Variant data;
}


struct ConnectorBase(T) {
  T[] inputs, outputs;
}

class NodeType {
  string name;
  ConnectorBase!NodeTypeParameter connections;
}

class NodeValue {
  ConnectorBase!NodeTypeParameter connections;
}

class NodeConnection {
  ConnectorBase!size_t connections;
}

class Node {
private:
  ImVec2 origin, size;
  size_t node_id;
  string name;
  NodeType node_type;
  NodeValue node_value;
  NodeConnection connection_ids;
public:
  this ( NodeType node_type_ ) {
    node_type = node_type_;
  }
}
