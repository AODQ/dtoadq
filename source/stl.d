module stl;
public import std.stdio, std.algorithm, std.math, std.string,
              std.array, std.bitmanip, std.concurrency, std.typecons,
              std.random, std.conv, std.exception;
public static import file = std.file,
                     regex = std.regex,
                     json = std.json;

/// --- move these elsewhere later ---

bool Not_Swap_File ( string file ) {
  auto rgx = stl.regex.ctRegex!r"\.#";
  return !stl.regex.matchFirst(file, rgx);
}


/** Returns if the file extension exists */
string RMatching_File_Extensions ( T... ) ( string str, T extensions ) {
  foreach ( ext; extensions ) {
    if ( str.length <= ext.length      ) continue;
    if ( str[$-ext.length .. $] == ext ) return ext;
  }
  return "";
}

string Truncate_Directory ( string filename ) {
  auto index = filename.lastIndexOf('/');
  if ( index == -1 ) return filename;
  return filename[index+1 .. $];
}

string Truncate_Extension ( string filename ) {
  auto index = filename.indexOf('.');
  if ( index == -1 ) return filename;
  return filename[0 .. index];
}

string Truncate_DirExt ( string filename ) {
  return filename.Truncate_Directory.Truncate_Extension;
}


import std.traits, std.meta;

/// returns all members of T
template AllMembers(alias T) {
  private template MemberFilter(string name) {
    mixin(`alias field = %s.%s;`.format(fullyQualifiedName!T, name));
    enum MemberFilter = !is(field) && field;
  }
  alias AllMembers = Filter!(MemberFilter, __traits(derivedMembers, T));
}

bool HasAttribute(alias T, alias member, string name)() {
  bool result = false; // To avoid "statement is not reachable"
                       //   ( Yes, really, wtf! )
  foreach ( attr; __traits(getAttributes, __traits(getMember, T, member)) ) {
    if ( attr.stringof == name ) result |= true;
  }
  return result;
}

// Returns string of all members of filtered type
string[] AllFilteredAttributes(alias T, string name)() {
  import functional;
  string[] results;
  foreach ( member; __traits(derivedMembers, T) ) {
    foreach ( attr; __traits(getAttributes, __traits(getMember, T, member)) ) {
      if ( attr.stringof == name ) {
        results ~= member;
        pragma(msg, "ADD: " ~ member);
      }
    }
  }
  return results;
}

/// RVal Reference (idiom written by Randy Schutt)
mixin template RValRef ( ) {
  static assert ( is(typeof(this) == struct) );

  ref const(typeof(this)) Ref() const pure nothrow return {
    return this;
  }
}
