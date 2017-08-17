module stl;
public import std.stdio, std.algorithm, std.math, std.string,
              std.array, std.bitmanip, std.concurrency, std.typecons,
              std.random, std.conv, std.exception;
public static import file = std.file,
                     regex = std.regex,
                     json = std.json;



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


/// Eponomyous template that returns all members of T
template AllMembers(alias T) {
  private template MemberFilter(string name) {
    mixin(`alias field = %s.%s;`.format(fullyQualifiedName!T, name));
    enum MemberFilter = !is(field) && field;
  }
  alias AllMembers = Filter!(MemberFilter, __traits(derivedMembers, T));
}


/// RVal Reference (idiom written by Randy Schutt)
mixin template RValRef ( ) {
  static assert ( is(typeof(this) == struct) );

  ref const(typeof(this)) Ref() const pure nothrow return {
    return this;
  }
}
