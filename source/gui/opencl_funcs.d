module gui.opencl_funcs;
import gui.tnodes : SubnodeDescription, SubnodeType, ParseNodeType, NodeType;
// module to generate node types and other meta funcs
private:
alias SD  = SubnodeDescription,
              ST  = SubnodeType,
              PNT = ParseNodeType,
              NT  = NodeType;

ST STV  ( string name ) { return ST (name, SD.Varying); }
ST STF3 ( string name ) { return ST (name, SD.Float3);  }
ST STF2 ( string name ) { return ST (name, SD.Float2);  }
ST STF  ( string name ) { return ST (name, SD.Float );  }
ST STC  ( string name ) { return ST (name, SD.Colour);  }
ST STI  ( string name ) { return ST (name, SD.Int   );  }
ST STS  ( string name ) { return ST (name, SD.String);  }

alias PFn = PNT.Function, POp = PNT.Operation, PCo = PNT.Constant;

NT NNT(string n, ST[] i, ST[] o, SD s ) { return new NT(n, i, o, s); }
NT NNT(string n, ST[] i, ST[] o, PNT p) { return new NT(n, i, o, p); }
NT Map_Func (T...)( string name, T args ) {
  ST[] stargs; foreach ( n; args ) stargs ~= ST (n);
  return new NT(name, [ST("Origin")] ~ stargs, [ST("Out")], PFn);
}

struct NodeCategory {
  string label;
  NodeType[] nodes;
}

static NodeCategory[] node_types = [
  NodeCategory("Data Types", [
    NNT("Colour", [], [STC(" "  )], SD.Colour ),
    NNT("Float3", [], [STF3(" " )], SD.Float3 ),
    NNT("Float2", [], [STF2(" " )], SD.Float2 ),
    NNT("Float",  [], [STF(" "  )], SD.Float  ),
    NNT("Int",    [], [STI(" "  )], SD.Int    ),
    NNT("String", [], [STS(" "  )], SD.String ),
  ]),
  NodeCategory("Variables", [
    NNT("Origin", [], [STF3(" " )], PCo ),
    NNT("Time"  , [], [STF(" "  )], PCo ),
  ]),
  NodeCategory("Mathematics", [
    NNT("+",   [STV("In0"  ), STV("In1"   )], [STV("Out" )], POp ),
    NNT("*",   [STV("In0"  ), STV("In1"   )], [STV("Out" )], POp ),
    NNT("-",   [STV("In0"  ), STV("In1"   )], [STV("Out" )], POp ),
    NNT("/",   [STV("In0"  ), STV("In1"   )], [STV("Out" )], POp ),
    NNT("cos", [STV("In" )],                  [STV("Out" )], PFn ),
    NNT("sin", [STV("In" )],                  [STV("Out" )], PFn ),
    NNT("tan", [STV("In" )],                  [STV("Out" )], PFn ),
  ]),
  NodeCategory("Map Basic", [
    NNT("Render_Map", [STF3("Function"), STI("Material")], [], PFn),
    Map_Func("SD_Sphere", "Radius"),
    Map_Func("SD_Box", "Bounds"),
    Map_Func("SD_Cross", "dist"),
  ]),
  NodeCategory("Map Mod", [
    Map_Func("OP_Union",     "Origin"),
    Map_Func("OP_Subtract",  "Origin"),
    Map_Func("OP_Intersect", "Origin"),
    Map_Func("OP_Repeat",    "Modulo"),
  ]),
];


// -- generator functions --
public auto RNode_Type_Map ( ) {
  import functional;
  NodeType[string] node_types_map;
  node_types.each!(cat => cat.nodes.each!(n => node_types_map[n.RName] = n));
  return node_types_map;
}

public void Apply_Function ( void delegate(string, string[]) fn ) {
  import functional;
  node_types.each!(category => fn(category.label,
                                  category.nodes.map!(n => n.RName).array));
}
