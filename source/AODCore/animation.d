module AODCore.animation;
import AODCore.aod;
import JSON = std.json;
import AODCore.aod;

/**
  An animation that stores an array of sheetrects
*/
class Animation {
public:
  /** Type of animation playback to occur */
  enum Type {
    /** 0 .. $ */
    Linear,
    /** $ .. 0 */
    Reverse,
    /** 0 .. $ .. 0 */
    Zigzag,
    /** */
    Nil
  }

  bool flip_x, flip_y;

  /** Animation textures in the order they should be played */
  SheetRect[] textures;
  /** frames that should occur per every texture iteration */
  int frames_per_texture;
  /** */
  Type type;
  /** Must be at least two textures */
  this(Type _type, SheetRect[] _textures, int _frames_per_texture) {
      flip_x = flip_y = false;
    type = _type;
    textures = _textures.dup;
    textures = _textures;
    frames_per_texture = _frames_per_texture;
  }
  /** Parses JSON object
    frames = "X, Y"
    order  = "Linear|Reverse|Zigzag" (if not supplied, Nil)
    Example: {"frames": "0..5, 1", "order": "linear"}
    Params:
      data   = data to to use to construct
      sc     = sheet container
      width  = sprite width
      height = sprite height
  **/
  this ( JSON.JSONValue data, SheetContainer sc, int width, int height ) {
    import std.string : strip;
    string sframes = data["frames"].str.strip;
    int sx, sy, ex = -1, ey = -1;
    { // extract animation indices
      int* digptr = &sx;
      foreach ( c; sframes ) {
        if ( c == ' ' ) continue;
        if ( c == '-' ) {
          if ( digptr == &sx ) {
            ex = 0;
            digptr = &ex;
          }
          if ( digptr == &sy ) {
            ey = 0;
            digptr = &ey;
          }
        }
        if ( c == ',' ) {
          if ( ex == -1 ) {
            ex = sx;
          }
          digptr = &sy;
        }
        if ( c >= '0' && c <= '9' ) {
          int digit = cast(int)(c - '0');
          *digptr = (*digptr * 10)+digit;
        }
      }
      if ( ey == -1 ) { ey = sy; }
      for ( int ix = sx; ix < ex+1; ++ ix ) {
      for ( int iy = sy; iy < ey+1; ++ iy ) {
        auto sr = SheetRect(sc, ix, iy, width, height);
        textures ~= sr;
      }}
    }

    { // extract type
      type = Type.Nil;
      if ( const(JSON.JSONValue)* order = ("order" in data) ) {
        import std.string : toLower;
        switch ( order.str.toLower ) {
          default: assert(0);
          case "linear"  : type = Type.Linear ; break;
          case "reverse" : type = Type.Reverse; break;
          case "zigzag"  : type = Type.Zigzag ; break;
        }
      }
    }

    { // extract frames per texture
      frames_per_texture = -1;
      if ( const(JSON.JSONValue)* fps = ("framerate" in data) ) {
        frames_per_texture = cast(int)To_MS(fps.integer);
      }
    }
  }
}

/**
  Keeps track of animation and time left along with index. Make sure to call
  Update_Index.
*/
struct Animation_Player {
public:
  /* @disable this(); */
  /** See Set */
  this(Animation _animation)
  { Set(_animation, true); }
  /** Animation to 'play' */
  Animation animation;
  /** Current texture index of the animation */
  size_t index;
  /** Frames left before index is incremented */
  int frames_left;
  /** The current direction of playback, 1 = positive.
      Only applies to Zigzag animation type */
  bool direction;
  /** Indicates whether the animation has finished playing, if another update
        is called after this is set, it is reset to false */
  bool done;
  /** Updates the animation player, should be called once every frame
Return:
    Current index
   */
  int Update() {
    done = false;
    if ( animation.textures.length <= 1 ) return 0;
    if ( -- frames_left <= 0 ) {
      frames_left = animation.frames_per_texture;
      switch ( animation.type ) {
        default: break;
        case Animation.Type.Linear:
          if ( ++ index >= animation.textures.length ) {
            index = 0;
            done = true;
          }
        break;
        case Animation.Type.Zigzag:
          if ( (direction && ++ index >= animation.textures.length) ) {
            index = animation.textures.length - 2;
            direction ^= 1;
            break;
          }
          if ( (!direction && -- index <  0) ) {
            index = 1;
            direction ^= 1;
          }
        break;
        case Animation.Type.Reverse:
          if ( -- index < 0 ) {
            done = true;
            index = animation.textures.length - 1;
          }
        break;
      }
    }
    return cast(int)index;
  }
  /** Returns current texture based off the index */
  SheetRect R_Current_Frame() {
    if ( animation.textures.length != 0 )
      return animation.textures[index];
    return SheetRect();
  }
  /** resets a new animation, identical to Set(animation, true) */
  void Reset(Animation _animation) {
    Set(_animation, true);
  }
  /** Sets a new animation
Params:
  _animation    = new animation
  force_reset   = if animation already playing, will force variables to reset
  override_type = overrides default animation type
    */
  void Set(Animation _animation, bool force_reset = false) {
    if ( animation !is _animation || force_reset ) {
      done = false;
      animation = _animation;
      frames_left = animation.frames_per_texture;
      direction = true;
      index = 0;
      if ( animation.type == Animation.Type.Reverse ) {
        index = animation.textures.length - 1;
        direction = false;
      }
    }
  }
  /** Returns if animation is null (no animation to play) */
  bool Valid ( ) { return animation !is null; }
}
