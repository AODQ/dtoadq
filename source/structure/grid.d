module structure.grid;
// import globals, ray;

// class Grid {
//   Voxel[][][] voxels;
//   size_t map_dim, voxel_dim;
//   size_t voxel_length;
// public:
//   this ( size_t voxel_length_, size_t voxel_dim_ ) {
//     voxel_length = voxel_length_;
//     voxels.length = voxel_length;
//     foreach ( y; voxels ) {
//       y.length = voxel_length;
//       foreach ( z; y )
//         z.length = voxel_length;
//     }
//     map_dim = cast(size_t)(voxel_dim * (voxel_length/2.0f));
//   }
//   void Insert ( Voxel voxel ) {
//     size_t x, y, z;
//     x = cast(size_t)voxel.origin.x;
//     y = cast(size_t)voxel.origin.y;
//     z = cast(size_t)voxel.origin.z;
//     // assert(voxels[x][y][z] is null);
//     voxels[x][y][z] = voxel;
//   }

//   auto RVoxel ( size_t x, size_t y, size_t z ) {return voxels[x][y][z];}
//   auto RIndices ( gln.vec3 position ) {
//     auto mapdim_vec = gln.vec3(map_dim, map_dim, map_dim);
//     return (position + mapdim_vec)/voxel_length;
//   }
//   void RIndices ( gln.vec3 position, out size_t x, out size_t y, out size_t z ){
//     auto ind = RIndices(position);
//     x = cast(size_t)ind.x;
//     y = cast(size_t)ind.y;
//     z = cast(size_t)ind.z;
//   }
//   auto RVoxel ( gln.vec3 position ) {
//     size_t x, y, z;
//     RIndices(position, x, y, z);
//     return voxels[x][y][z];
//   }
// }

// auto Construct_Grid ( 

