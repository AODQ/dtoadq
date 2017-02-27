module image;
static import img = imageformats;


float[] Read_Image(string filename) {
  import functional;
  ubyte[] data = img.read_tga(filename).pixels;
  float[] fdata = data.map!(n => (n/255.0f)).array;
  return fdata;
}
