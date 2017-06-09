#!/usr/bin/env rdmd
import std.stdio;
import std.string;

class Dir {
  string dir_name;
  Dir[string] dirs;
  string[] files;

  this ( string dir_name_ ) { dir_name = dir_name_; }

  void Add ( bool is_file, string fname, string name = "" ) {
    if ( name == "" )
      return this.Add(is_file, fname, fname[fname.indexOf("/")+1..$]);
    auto ind = name.indexOf("/");
    if ( ind == -1 ) {
      if ( is_file ) files ~= name;
      else dirs[name] = new Dir(fname);
    } else {
      string dir = name[0..ind];
      writeln("DIR: ", dir);
      writeln("NAM: ", name);
      writeln("FIL: ", name[ind+1..$]);
      writeln("++");
      this.Print();
      writeln("++");
      writeln("+");
      dirs[dir].Add(is_file, fname, name[ind+1..$]);
    }
  }
  void Print ( string spaces = "" ) {
    writeln(spaces, "D ", dir_name);
    spaces ~= "  ";
    foreach ( f; files )
      writeln(spaces, "F ", f);
    foreach ( d; dirs )
      d.Print(spaces);
  }
}

void main (  ) {
  static import Fil = std.file;
  import std.algorithm, std.array, std.string;
  string dir = "projects/";
  Dir base_dir = new Dir(dir);
  foreach ( fil; Fil.dirEntries(dir, Fil.SpanMode.breadth) ) {
    base_dir.Add(fil.isFile, fil.name);
  }

  writeln("----");
  base_dir.Print();
}

