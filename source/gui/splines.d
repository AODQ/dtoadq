module gui.splines;
import derelict.imgui.imgui;

void Draw_Hermite ( ImDrawList* draw_list, ImVec2 p1, ImVec2 p2 ) {
  ImVec2 v1 = ImVec2(80.0f, 0.0f), v2 = v1;
  size_t steps = 20;
  foreach ( step; 0 .. steps ) {
    float t = cast(float)step/cast(float)steps, t2 = t*t, t3 = t*t*t;
    float h1 =  2*t3 - 3*t2 + 1.0f, h2 = -2*t3 + 3*t2,
          h3 =    t3 - 2*t2 + t,    h4 =    t3 -   t2;
    ImDrawList_PathLineTo(draw_list,
      ImVec2(h1*p1.x + h2*p2.x + h3*v1.x + h4*v2.x,
             h1*p1.y + h2*p2.y + h3*v1.x + h4*v2.y));
  }
  ImDrawList_PathStroke(draw_list, ImColor(0.3f, 0.2f, 0.6f), false, 3.0f);
}
