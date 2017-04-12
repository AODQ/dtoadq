module softwarerenderer.image;
import globals;
import derelict.opengl3.gl3;

class Image {
  float[] pixels;
  size_t width, height;
  GLuint texture;
public:
  this ( size_t width_, size_t height_ ) {
    width = width_; height = height_;
    pixels.length = height*height*3;
    import std.random;
    foreach ( ref p; pixels ) p = uniform(0.0f, 1.0f);
    glGenTextures(1, &texture);
  }

  void Apply ( size_t x, size_t y, gln.vec3 pixel ) {
    pixels[(y*width + x)*3 + 0] = pixel.r;
    pixels[(y*width + x)*3 + 1] = pixel.g;
    pixels[(y*width + x)*3 + 2] = pixel.b;
  }

  AOD.SheetContainer To_OGL_Sprite ( ) {
    glGetError();
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, cast(int)width, cast(int)height, 0,
        GL_RGB, GL_FLOAT, cast(void*)pixels.ptr
    );
    glBindTexture(GL_TEXTURE_2D, 0);
    texture.writeln;
    return AOD.SheetContainer(texture, cast(int)width, cast(int)height);
  }
}
