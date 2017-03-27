module cloctree;
import globals;
import opencl;
/**
  In order to fit the octree on the GPU, I have to use a pool, so the
  data must be linear.
*/

/**
  Either the node of an octree or the octree itself. To get the voxel/children
  use an OctreeData, preferrably with the RNode and RVoxel functions. Note that
  the lack of an ID - the equivalence of a null - is -1.

  Also contains information on the origin of the node as well as half its size
    so to find its minimum X point, origin[0] - half_size[0]
*/
struct CLOctreeNode {
  cl_int8 child_id;
  cl_float3 origin;
  cl_float3 half_size;
  cl_int voxel_id;
  cl_ulong buffer, buffer2, buffer3;
}

/** */
auto New_CLOctreeNode ( inout float[3] origin_, inout float[3] half_size_ ) {
  CLOctreeNode node = {
    origin    : To_CLFloat3(origin_),
    half_size : To_CLFloat3(half_size_),
    child_id  : To_CLInt8([-1, -1, -1, -1, -1, -1, -1, -1]), voxel_id : -1,
  };
  return node;
}


/**
  Contains BRDF, BTDF, size and position data
*/
struct CLVoxel {
  cl_float3 position;
  cl_float3 colour;
}

/** */
auto New_CLVoxel ( float[3] position_ ) {
  import std.random;
  CLVoxel voxel = {
    position : To_CLFloat3(position_),
    colour : To_CLFloat3([uniform(0.5f, 1.0f), uniform(0.0f, 1.0f),
                          uniform(0.5f, 1.0f)])
  };
  return voxel;
}

/**
  Contains data of the nodes/voxels
*/
struct OctreeData {
  CLOctreeNode[] node_pool;
  CLVoxel[]      voxel_pool;
}

/**
  Returns a copy of the node
*/
CLOctreeNode RNode(inout OctreeData data, int node_id) {
  return *data.RNodePtr(node_id);
}

/**
  Returns a pointer to the node
*/
auto RNodePtr(inout OctreeData data, int node_id) in {
  assert(node_id >= 0 && node_id < data.node_pool.length,
         "node id out of range: " ~ node_id.to!string);
}   body   {
  return &data.node_pool[cast(size_t)node_id];
}

/**
  returns a copy of the voxel
*/
CLVoxel RVoxel(inout OctreeData data, int voxel_id) {
  return *data.RVoxelPtr(voxel_id);
}
/**
  returns a pointer to the voxel
*/
auto RVoxelPtr(inout OctreeData data, int voxel_id) in {
  assert(voxel_id >= 0 && voxel_id < data.voxel_pool.length,
         "voxel id out of range: " ~ voxel_id.to!string);
}   body   {
  return &data.voxel_pool[cast(size_t)voxel_id];
}

/**
  Constructs an octree from an origin and original dimensions (note its half
    of the preferred dimensions), requires you fill out the voxel data
    beforehand.
*/
auto Construct_CLOctree ( float[3] origin, float[3] half_size,
                          inout CLVoxel[] data ) {
  OctreeData octree_data = {
    node_pool : [ New_CLOctreeNode(origin, half_size) ],
    voxel_pool : data.dup
  };
  foreach ( i; 0 .. data.length ) octree_data.Insert(cast(int)i);
  return octree_data;
}

/**
  Returns the octant that the node lies in, more specifically, it returns
  the bitmask according to this diagram:
        ID : 0 1 2 3 4 5 6 7
  point > x: f f f f T T T T
  point > y: f f T T f f T T
  point > z: f T f T f T f T
*/
ubyte ROctant_Mask ( inout CLOctreeNode node, inout cl_float4 point ) {
  import functional;
  ubyte oct = 0;
  foreach ( i; 0 .. 3 )
    oct |= 4/(1+i)*(node.origin[i] < point[i]);
  return oct;
}


/**
  Returns if the node is a leaf, that is, it has no children.
*/
bool Is_Leaf ( inout CLOctreeNode node ) {
  return node.child_id[0] == -1;
}

/**
  Returns true if node is empty (has no children nor data)
*/
bool Is_Empty (inout CLOctreeNode node) {
  return node.Is_Leaf && node.voxel_id == -1;
}

/**
  Inserts the node into the tree, given the voxel exists already
*/
void Insert ( ref OctreeData data, int voxel_id, int node_id = 0 ) in {
  assert(voxel_id >= 0 && voxel_id < data.voxel_pool.length,
         "voxel id out of range: " ~ voxel_id.to!string);
  assert(node_id >= 0 && node_id < data.node_pool.length,
         "node id out of range: " ~ node_id.to!string);
}   body   {
  auto node = &data.node_pool[cast(size_t)node_id];
  if ( (*node).Is_Leaf ) {
    if ( node.voxel_id == -1 ) {
      // if it's a leaf and the voxel ID is not set, we can set it for now
      node.voxel_id = voxel_id;
      return;
    } else {
      // Not enough room for two IDs, have to create a new set of octants
      // and insert the IDs into the octants
      int old_voxel_id = node.voxel_id;
      node.voxel_id = -1;
      int[8] child_id;

      foreach ( i; 0 .. 8 ) {
        import functional;
        float[3]
             new_origin = [node.origin[0], node.origin[1], node.origin[2]],
             new_dim    = node.half_size.array.map!"a*0.5f".array[0..3]
                              .to!(float[3]);
        foreach ( p; 0 .. 3 )
          new_origin[p] += node.half_size[p]* (i&(4/(1+p)) ? 0.5f : -0.5f);
        child_id[i] = cast(int)data.node_pool.length;
        data.node_pool ~= New_CLOctreeNode(new_origin, new_dim);
      }

      // I get some weird problem if I just use node.child_id[i] in
      // the above for loop; only the first element persists in the node pool
      node.child_id = child_id.dup;

      auto old_index = (*node).ROctant_Mask(data.RVoxel(old_voxel_id).position),
           new_index = (*node).ROctant_Mask(data.RVoxel(    voxel_id).position);
      Insert(data, old_voxel_id, node.child_id[old_index]);
      Insert(data,     voxel_id, node.child_id[new_index]);
    }
  } else {
    // Just recursively insert the node into the corresponding child
    // until we hit a leaf node
    auto index = (*node).ROctant_Mask(data.RVoxel(voxel_id).position);
    assert(node !is &data.node_pool[node.child_id[index]],
           "recursive node insertion - degenerate tree formed");
    Insert(data, voxel_id, node.child_id[index]);
  }
}

/**


  .------.
  |\   7 |\
  | '------'
  |1|  5 |3|
  '-|----' | // front face is 2
   \|  6  \|
    '------'

        ID : 0 1 2 3 4 5 6 7
  point > x: f f f f T T T T
  point > y: f f T T f f T T
  point > z: f T f T f T f T
*/

struct RopeInfo {
  int id;
  ubyte mask;
}

/** Calculate face

  I calculate neighbouring nodes in the parent octree for trivial
  cases:

  .--.--.
  | 2| 6|
  [--+--]
 1| 0| 4|
  '--'--'

  thus for 0, at faces 7, 3, 5 I know the siblings 2, 4, 1 respectively

  if it is not known, I go up a parent and calculate. I can't simply take
  my neighbour as one might not exist; this is a sparse octree. So taking
  my "uncle" is close enough, at least for now

  if I run out of parents, then the node is an outside boundary node with
  a face facing outwards
*/
int Calculate_Face ( ref OctreeData data, RopeInfo[] info, ubyte face, ubyte mask ) {
  import std.algorithm.mutation : reverse;
  if ( info.length == 1 ) return -1;
  auto parent = parents[$-1];
  switch ( mask ) {
    case 0:
      switch ( face ) {
        case 7: return data.RNode(parent.id).child_id[2];
        case 3: return data.RNode(parent.id).child_id[4];
        case 5: return data.RNode(parent.id).child_id[1];
        default: return data.Calculate_Face(info[1 .. $-1], 0);
      }
    break;
  }
  return -1;
}

/** Calculate rope */
void Calculate_Ropes ( ref OctreeData data, RopeInfo[] info = [],
                       int node_id = 0, ubyte node_mask = 0 ) {
  auto node = data.RNode(node_id);
  if ( node.Is_Leaf && node.voxel_id != -1 ) {
    foreach ( i; 0 .. 6 ) {
      auto ui = cast(ubyte)i;
      data.RNode(node_id).child_id[i+1] = data.Calculate_Face(parents, ui);
    }
  } else {
    info ~= RopeInfo(node_id, node_mask);
    foreach ( i; 0 .. 8 )
      Calculate_Ropes(data, parents, i, node.child_id[i]);
  }
}


/**
  Counts the amount of nodes in the tree, mostly to check for degeneracy.
  It's much faster to just use data.node_pool.length
*/
int Count_Nodes ( inout OctreeData data, int node_id = 0 ) {
  auto node = data.RNode(node_id);
  if ( node.Is_Leaf ) {
    return cast(int)(node.voxel_id != -1);
  } else {
    int 서 = 0;
    foreach ( i; 0 .. 8 )
      서 += Count_Nodes(data, node.child_id[i]);
    return 서;
  }
}

/**
  Sets min/max as the bounds of the node
*/
void RBounds ( inout OctreeData data, int node_id,
               out float[3] min, out float[3] max ) {
  return data.RNode(node_id).RBounds(min, max);
}

/**
  Returns the bounds of the node as [min, max]
*/
float[3][2] RBounds ( inout CLOctreeNode node ) {
  float[3] min, max;
  RBounds(node, min, max);
  return [min, max];
}

/**
  Sets min/max as the bounds of the node
*/
void RBounds ( inout CLOctreeNode node, out float[3] min, out float[3] max ) {
  foreach ( i; 0 .. 3 ) {
    max[i] = node.origin[i] + node.half_size[i];
    min[i] = node.origin[i] - node.half_size[i];
  }
}

/**
  Returns as list of voxel ids for all voxels within the box
  described by min .. max
*/
int[] RVoxels_Inside_Box ( inout OctreeData data, float[3] min, float[3] max,
                           int node_id = 0) {
  auto node = data.RNode(node_id);
  if ( node.Is_Leaf ) {
    if ( node.voxel_id != -1 ) {
      float[] p = data.RVoxel(node.voxel_id).position;
      if ( p[0] > max[0] || p[1] > max[1] || p[2] > max[2] ) return [];
      if ( p[0] < min[0] || p[1] < min[1] || p[2] < min[2] ) return [];
      return [ node.voxel_id ];
    }
    return [];
  } else {
    int[] results;
    foreach ( i; 0 .. 8 ) {
      float[3] cmax, cmin;
      RBounds(data, node_id, cmin, cmax);

      if ( cmax[0] < min[0] || cmax[1] < min[1] || cmax[2] < min[2]) continue;
      if ( cmin[0] > max[0] || cmin[1] > max[1] || cmin[2] > max[2]) continue;

      results ~= data.RVoxels_Inside_Box(min, max, node.child_id[i]);
    }
    return results;
  }
}

unittest {
  import std.stdio;
  float Rand() {
    import std.random;
    return uniform(-1.0f, 1.0f);
  }

  float[3] Rand_Vec() {
    return [Rand(), Rand(), Rand()];
  }

  bool Point_In_Box ( float[3] point, float[3] min, float[3] max ) {
    return point[0] >= min[0] && point[1] >= min[1] && point[2] >= min[2] &&
           point[0] <= max[0] && point[1] <= max[1] && point[2] <= max[2];
  }

  float[3][] points;
  stdout.flush();

  const size_t amt_points = 1_000_000;
  foreach ( i; 0 .. amt_points ) {
    points ~= Rand_Vec();
  }

  CLVoxel[] voxels;
  foreach ( ref p; points ) {
    voxels ~= New_CLVoxel(p);
  }

  float[3] tree_origin = [0.0f, 0.0f, 0.0f],
           tree_halef_siz = [1.0f, 1.0f, 1.0f];
  import functional;
  auto tree = Construct_CLOctree(tree_origin, tree_halef_siz, voxels);
  // query box
  float[3] qmin = [-0.05f, -0.05f, -0.05f],
           qmax = [ 0.05f,  0.05f,  0.05f];

  // -- asserting insertion and bounds testing
  import std.datetime;
  {
    size_t Bruteforce_Test ( ) {
      int count;
      foreach ( i; 0 .. amt_points )
        count += cast(int)(Point_In_Box(points[i], qmin, qmax));
      return count;
    }

    size_t Octree_Test ( ) {
      return tree.RVoxels_Inside_Box(qmin, qmax).length;
    }

    auto result = benchmark!(Bruteforce_Test, Octree_Test)(25);
    auto bf_result = result[0].msecs,
         tt_result = result[1].msecs;
    writeln("Octree Unittest Results, ", amt_points, " points");
    writeln("Query box: ", qmin, " - ", qmax);
    writeln("Bruteforce time: ", bf_result, " milliseconds");
    writeln("Octree     time: ", tt_result, " milliseconds");
    assert(Bruteforce_Test == Octree_Test);
  }
}
