module AODCore.serializer;
import std.traits;
import std.typetuple;
import std.typecons;
import std.meta;
import std.string : format;
import AOD = AODCore.aod;
import std.stdio : writeln;

/**
  Provides basic serialization and debugging tools
  See example unit test classes to see how it works.

  Current limitations:
    - serialization does not work on members that are classes or structs
  Weird limitations (never will be fixed):
    - serialization will not work on serialized members that are non-serialized
      classes
    - function overloading will only take length into consideration, not
      parameter. (it will work from derived class to the parent)
**/

/// mixin to enable serialization/deserialization on a class that inherits
/// another serializeable
mixin template SerializeClassDerived(alias T) {
  override void Serialize(ref std.json.JSONValue value) {
    import std.conv : to;
    mixin(AOD.Serializer.SerializeMembersMixin!T);
    super.Serialize(value);
  }
  override void Deserialize(std.json.JSONValue value) {
    import std.conv : to;
    mixin(AOD.Serializer.DeserializeMembersMixin!T);
    super.Deserialize(value);
    Deserialize_Extra();
  }
  override size_t R_Sizeof ( ) {
    return __traits(classInstanceSize, T);
  }
}

/// mixin to enable serialization/deserialization on a base class (no
/// serialization on any inherited class)
mixin template SerializeClassBase(alias T) {
  void Serialize(ref std.json.JSONValue value) {
    import std.conv : to;
    static import AODCore.serializer;
    mixin(AOD.Serializer.SerializeMembersMixin!T);
  }
  void Deserialize(std.json.JSONValue value) {
    import std.conv : to;
    mixin(AOD.Serializer.DeserializeMembersMixin!T);
    Deserialize_Extra();
  }
  size_t R_Sizeof ( ) {
    return __traits(classInstanceSize, T);
  }
  void Deserialize_Extra() {}
}

/// FIXME !!! doesn't work ?!
/// Supposed to be identical to DebugFunctionClass*, yet doesn't work
mixin template DebugFunctionClassDerived__BROKEN(string T) {
  override string[] R_Debug_Function_Names_Collector() {
    mixin(`return [AllDebugFuncs!%s]`.format(T));
  }
}

/// mixin to enable debug functions on a base class
string DebugFunctionClassMixinBase(string T)() {
  return `
    string[] R_Debug_Function_Names() {
      return [AllDebugFuncs!%s];
    }
  `.format(T);
}

/// mixin to enable both serialize and debug functions on a derived class
string SerializeDebugMixinDerived(string T)() {
  return format(q{
    mixin AOD.Serializer.SerializeClassDerived!%s;
    mixin(`
      override string[] R_Debug_Function_Names() {
        return [AOD.Serializer.AllDebugFuncs!%s] ~ super.R_Debug_Function_Names;
      }
    `);
    override bool Debug_Function_Call(string fn, string[] args) {
      import std.conv : to;
      mixin(AOD.Serializer.CallFunctions!%s);
      return super.Debug_Function_Call(fn, args);
    }
  }, T, T, T);
}

/// mixin to enable both serialize and debug functions on a base class
string SerializeDebugMixinBase(string T)() {
  return `
    mixin AOD.Serializer.SerializeClassBase!%s;
    string[] R_Debug_Function_Names() {
      return [AOD.Serializer.AllDebugFuncs!%s];
    }
    bool Debug_Function_Call(string fn, string[] args) {
      mixin(AOD.Serializer.CallFunctions!%s);
      return false;
    }
  `.format(T, T, T);
}

// -----------------------------------------------------------------------------
// ------------------------------ utilities ------------------------------------
// -----------------------------------------------------------------------------

// Returns all members that are serializeable (have attrib tag), does not check
// and inherited classes
template AllSerializeableMembers(alias T) {
  private template SerializeableMember(string attrib) {
    enum SerializeableMember = attrib == "serialize";
  }
  private template MemberFilter(string name) {
    mixin(`alias field = %s.%s;`.format(fullyQualifiedName!T, name));
    enum MemberFilter = !is(field) && !isCallable!field &&
            anySatisfy!(SerializeableMember, __traits(getAttributes, field));
  }
  alias AllSerializeableMembers =
               Filter!(MemberFilter, __traits(derivedMembers, T));
}

// AllSerializeableMembers unittest
unittest {
  import std.algorithm : sort;
  assert([AllSerializeableMembers!UnitTestA] == []);
  assert([AllSerializeableMembers!UnitTestB].sort() == ["x", "y", "c"].sort());
  assert([AllSerializeableMembers!UnitTestC].sort() == ["f", "z"].sort());
  assert([AllSerializeableMembers!UnitTestD] == []);
}

// Returns the types of members generated from AllSerializeableMembers
template AllSerializeableMemberTypes(alias T) {
  private template MemberType(string name) {
    mixin(`alias field = %s.%s;`.format(fullyQualifiedName!T, name));
    enum MemberType = typeof(field).stringof;
  }
  alias AllSerializeableMemberTypes = staticMap!(MemberType,
                                                 AllSerializeableMembers!T);
}

// AllSerializeableMemberTypes unittest
unittest {
  import std.algorithm : sort;
  assert([AllSerializeableMemberTypes!UnitTestA] == []);
  assert([AllSerializeableMemberTypes!UnitTestB].sort() ==
                          ["int", "int", "bool"].sort());
  assert([AllSerializeableMemberTypes!UnitTestC].sort() ==
                              ["string", "uint"].sort());
  assert([AllSerializeableMemberTypes!UnitTestD] == []);
}

// Returns all debuggable functions as an AliasSeq for this derived
template AllDebugFuncs(alias T) {
  private template DebuggeableFunc(string attrib) {
    enum DebuggeableFunc = attrib == "debug";
  }
  private template FunctionFilter(string name) {
    mixin(`alias field = %s.%s;`.format(fullyQualifiedName!T, name));
    enum FunctionFilter = !is(field) && isCallable!field &&
              anySatisfy!(DebuggeableFunc, __traits(getAttributes, field));
  }
  alias AllDebugFuncs = Filter!(FunctionFilter, __traits(derivedMembers, T));
}

/// AllDebugFuncs unittest
unittest {
  import std.algorithm : sort;
  assert([AllDebugFuncs!UnitTestA].length == 0);
  assert([AllDebugFuncs!UnitTestB].sort() == ["thing2", "thing"].sort());
  assert([AllDebugFuncs!UnitTestC].length == 0);
  assert([AllDebugFuncs!UnitTestD].sort() == ["thing", "thing3"].sort());
}


// returns all function parameters as an AliasSeq
template AllFunctionParameters(alias T) {
  private template FunctionType(string name) {
    mixin(`alias field = %s.%s;`.format(fullyQualifiedName!T, name));
    alias FunctionType = AODCore.util.AliasSeqToStringArray!(Parameters!field);
  }
  alias AllFunctionParameters = staticMap!(FunctionType, AllDebugFuncs!T);
}

/// AllFunctionParameters unittest
unittest {
  import std.algorithm : sort;
  assert([AllFunctionParameters!UnitTestA].length == 0);
  assert([AllFunctionParameters!UnitTestB].sort() ==
                [["int"], ["int", "string", "ulong"]].sort());
  assert([AllFunctionParameters!UnitTestC].length == 0);
  assert([AllFunctionParameters!UnitTestD].sort() ==
                [["bool", "string"], []].sort());
}

/// ----------------------------------------------------------------------------
/// ------------------------- magic tricks -------------------------------------
/// ----------------------------------------------------------------------------

private:
struct MemberInfo { string name, type; }
struct FuncInfo { string name; string[] param_types; }

MemberInfo[] R_MemberInfo(T)() {
  MemberInfo[] members;
  string[] member_name   = [AllSerializeableMembers!T         ],
           member_types  = [AllSerializeableMemberTypes!T     ];
  foreach ( i; 0 .. member_name.length )
    members ~= MemberInfo(member_name[i], member_types[i]);
  return members;
}

// R_MemberInfo unittest
unittest {
  import std.algorithm : sort;
  assert(R_MemberInfo!UnitTestA == []);
  assert(R_MemberInfo!UnitTestB == [
    MemberInfo("x", "int"),
    MemberInfo("y", "int"),
    MemberInfo("c", "bool")
  ]);
  assert(R_MemberInfo!UnitTestC == [
    MemberInfo("f", "string"),
    MemberInfo("z", "uint")
  ]);
  // assert(R_MemberInfo!UnitTestD == [
  //   MemberInfo("a", "UnitTestA"),
  //   MemberInfo("b", "UnitTestB")
  // ]);
}

FuncInfo[] R_FuncInfo(T)() {
  FuncInfo[] funcs;
  string[]   func_name  = [AllDebugFuncs!T];
  string[][] func_param = [AllFunctionParameters!T];
  assert(func_name.length == func_param.length);
  foreach ( i; 0 .. func_name.length )
    funcs ~= FuncInfo(func_name[i], func_param[i]);
  return funcs;
}

// R_FuncInfo unittest
unittest {
  import std.algorithm : sort;
  assert(R_FuncInfo!UnitTestA == []);
  assert(R_FuncInfo!UnitTestB == [
    FuncInfo("thing2", ["int"]),
    FuncInfo("thing",  ["int", "string", "ulong"])
  ]);
  assert(R_FuncInfo!UnitTestC == []);
  assert(R_FuncInfo!UnitTestD == [
    FuncInfo("thing", ["bool", "string"]),
    FuncInfo("thing3", [])
  ]);
}

public string SerializeMembersMixin(T)() {
  string ret = "";
  foreach ( m; R_MemberInfo!T ) {
    ret ~= `value.object["%s"] = to!string(%s);`
           .format(m.name, m.name);
  }
  return ret;
}

public string DeserializeMembersMixin(T)() {
  string ret = "";
  foreach ( m; R_MemberInfo!T ) {
    ret ~= `
      if ( "%s" in value ) {
        try {
          %s = to!%s(value["%s"].str);
        } catch ( Exception e ) {
          writeln("invalid conversion");
        }
      }
    `.format(m.name, // if
             m.name, m.type, m.name // set value
     );
  }
  return ret;
}

/// Generates code so you can call "serialized" functions of type T
/// Requires the arguments of the function as an array of strings "args"
/// and the name of the function "fn" to be in scope
public string CallFunctions(T)() {
  auto funcs = R_FuncInfo!T;
  string Format_Params(string[] types) {
    string ret = "";
    foreach ( it; 0 .. types.length ) {
      ret ~= `to!(%s)(args[%s]),`.format(types[it], it);
    }
    if ( ret != "" ) ret = ret[0 .. $-1];
    return ret;
  }

  string ret = `
    import std.stdio : writeln;
    import std.conv : to;
    switch ( fn~to!string(args.length) ) {
      default: break; `;
  // the following foreach fills in the rest of the switch from the above ret
  foreach ( func; funcs ) {
    string fn_name = '"' ~ func.name ~ '"';
    ret ~= `
      case %s~"%s": // fn_name, fn_param_types.length
          if ( %s.length != args.length ) { // params
            writeln("parameter mismatch, needs: ", %s.length,  // params
                                   ", given: ", args.length);
            break;
          }
          try {
            mixin(%s ~ "(%s);"); // fn name, formatted params
            return true;
          } catch ( Exception e ) {
            writeln("invalid conversion");
          }
        break;
      `.format(fn_name, func.param_types.length,
               func.param_types, func.param_types,
               fn_name, Format_Params(func.param_types));
  }
  ret ~= "}";
  return ret;
}

// Test serialize, deserialize and call functions
unittest {
  UnitTestA A = new UnitTestA;
  UnitTestB B = new UnitTestB;
  UnitTestC C = new UnitTestC;
  UnitTestD D = new UnitTestD;
  // --- call functions
  writeln("1 A");
  assert(A.Debug_Function_Call("skipme", ["1"])            == false);
  writeln("1000 B");
  assert(B.Debug_Function_Call("thing2", ["1000"])         == true );
  writeln("invalid B");
  assert(B.Debug_Function_Call("thing2", ["inva"])         == false);
  assert(B.Debug_Function_Call("thing2", [])               == false);
  assert(B.Debug_Function_Call("thing2", ["true", "true"]) == false);
  assert(B.Debug_Function_Call("thing2", ["invalid"])      == false);
  assert(C.Debug_Function_Call("thing2", ["15"])           == true );
  assert(C.Debug_Function_Call("skipme", [])               == false);
  assert(D.Debug_Function_Call("thing",  ["true", "asdf"]) == true );
  assert(D.Debug_Function_Call("thing",  ["15", "f", "0"]) == true );
  assert(D.Debug_Function_Call("thing3", [])               == true );
  // --- deserialize
  // --- serialize
}

// -----------------------------------------------------------------------------
// ------------------------------ unit test classes ----------------------------
// -----------------------------------------------------------------------------
/// tests skipmes
private class UnitTestA {
public:
  this() {}
  int skip_me;
  void skipme() {}
  mixin(SerializeDebugMixinBase!"UnitTestA");
}
/// tests basic features
private class UnitTestB : UnitTestA {
public:
  this() {}
  @("serialize") {
    int x, y;
    bool c;
  }
  @("debug") void thing2(int) {}
             void skipme2 ( ) {}
  @("debug") void thing(int, string, ulong) {}

  mixin(SerializeDebugMixinDerived!"UnitTestB");
}
/// tests inheritance (can we still access B from C?)
private class UnitTestC : UnitTestB {
public:
  this() {}
  @("serialize") {
    string f;
    uint z;
  }
  mixin(SerializeDebugMixinDerived!"UnitTestC");
}
/// tests advanced features (overloading, empty functions, serializeable classes
private class UnitTestD : UnitTestC {
public:
  this() {}
  @("debug") void thing(bool, string) {}
  @("debug") void thing3() {}

  // TODO
  // @("serialize") {
  //   UnitTestA a;
  //   UnitTestB b;
  // }
  mixin(SerializeDebugMixinDerived!"UnitTestD");
}
