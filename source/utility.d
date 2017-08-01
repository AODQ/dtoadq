module utility;
import globals;
import std.string;

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
