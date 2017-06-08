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
