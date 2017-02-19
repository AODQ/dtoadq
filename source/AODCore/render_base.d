module AODCore.render_base;

import AODCore.aod;

/**
  The base of rendering, I suppose if you wanted to create your own entity,
  text, etc for some reason you could use this.
*/
class Render_Base {
public:
  this(ubyte _layer = 5) {
    layer = _layer;
  }
  abstract Vector R_Position(bool apply_static = false);
  /** Is the render_base visible? (is it rendered) */
  bool R_Visible() { return true;  }
  /** Returns layer */
  ubyte R_Layer()  { return layer; }
  /** Called whenever AOD adds this to the realm (no reason to have a
      Removed_From_Realm since doing so calls the destructor) */
  void Added_To_Realm() { }
  /** Does rendering to the screen */
  abstract void Render();
  /** Called once every frame. Meant to be overriden. */
  abstract void Update();
protected:
private:
  /** The layer (z-index) of which the object is located. Used only to
      determine which objects get rendered first */
  ubyte layer;
}
