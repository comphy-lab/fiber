/**
 * @brief Configures the camera view and rendering parameters for visualization.
 *
 * Parses and applies camera position, orientation, scaling, projection, image size, background color, and other rendering options to set up the viewing transformation for graphical output.
 */
else if (!strcmp (s, "view")) {
  float tx = 0.;
  float ty = 0.;
  float fov = 0.;
  float quat[4] = {0};
  float sx = 1.;
  float sy = 1.;
  float sz = 1.;
  unsigned width = 800;
  unsigned height = 800;
  unsigned samples = 4;
  float bg[3] = {0};
  float theta = 0.;
  float phi = 0.;
  float psi = 0.;
  bool relative = false;
  float tz = 0.;
  float near = 0.;
  float far = 0.;
  float res = 0.;
  char * camera = NULL;
  MapFunc map = NULL;
  int cache = 0;
  float p1x = 0.;
  float p1y = 0.;
  float p2x = 0.;
  float p2y = 0.;
  bview * view1 = NULL;
  Params params[] = {
    {"tx", pfloat, &tx},
    {"ty", pfloat, &ty},
    {"fov", pfloat, &fov},
    {"quat", pfloat, quat, 4},
    {"sx", pfloat, &sx},
    {"sy", pfloat, &sy},
    {"sz", pfloat, &sz},
    {"width", punsigned, &width},
    {"height", punsigned, &height},
    {"samples", punsigned, &samples},
    {"bg", pfloat, bg, 3},
    {"theta", pfloat, &theta},
    {"phi", pfloat, &phi},
    {"psi", pfloat, &psi},
    {"relative", pbool, &relative},
    {"tz", pfloat, &tz},
    {"near", pfloat, &near},
    {"far", pfloat, &far},
    {"res", pfloat, &res},
    {"camera", pstring, &camera},
    {"cache", pint, &cache},
    {"p1x", pfloat, &p1x},
    {"p1y", pfloat, &p1y},
    {"p2x", pfloat, &p2x},
    {"p2y", pfloat, &p2y},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  view (tx,ty,fov,quat,sx,sy,sz,width,height,samples,bg,theta,phi,psi,relative,tz,near,far,res,camera,map,cache,p1x,p1y,p2x,p2y,view1);
}
/**
 * @brief Handles the "draw_vof" command to render a volume-of-fluid (VOF) interface visualization.
 *
 * Parses parameters for VOF interface rendering, including component and surface names, edge display, fill mode, color settings, value ranges, colormap, face and line colors, line width, and expression mode. Invokes the draw_vof() function with the parsed arguments. Returns false if parameter parsing or rendering fails.
 */
else if (!strcmp (s, "draw_vof")) {
  char * c = NULL;
  char * s = NULL;
  bool edges = false;
  double larger = 0.;
  int filled = 0;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1.;
  bool expr = false;
  Params params[] = {
    {"c", pstring, &c},
    {"s", pstring, &s},
    {"edges", pbool, &edges},
    {"larger", pdouble, &larger},
    {"filled", pint, &filled},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!draw_vof (c,s,edges,larger,filled,color,min,max,spread,linear,map,fc,lc,lw,expr))
    return false;
}
else if (!strcmp (s, "isoline")) {
  char * phi = NULL;
  double val = 0.;
  int n = 1;
  bool edges = false;
  double larger = 0.;
  int filled = 0;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1.;
  bool expr = false;
  Params params[] = {
    {"phi", pstring, &phi},
    {"val", pdouble, &val},
    {"n", pint, &n},
    {"edges", pbool, &edges},
    {"larger", pdouble, &larger},
    {"filled", pint, &filled},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!isoline (phi,val,n,edges,larger,filled,color,min,max,spread,linear,map,fc,lc,lw,expr))
    return false;
}
else if (!strcmp (s, "cells")) {
  coord n = {0,0,1};
  double alpha = 0.;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!cells (n,alpha,lc,lw))
    return false;
}
/**
 * @brief Renders a vector field visualization with customizable appearance.
 *
 * Parses parameters for the vector field name, scaling factor, line color, line width, and level, then draws the vectors using these settings.
 *
 * @param u Name of the vector field to visualize.
 * @param scale Scaling factor for vector magnitudes.
 * @param lc RGB color array for vector lines.
 * @param lw Line width for vector rendering.
 * @param level Grid or refinement level to visualize; -1 for default.
 * @return false if parameter parsing or vector rendering fails, true otherwise.
 */
else if (!strcmp (s, "vectors")) {
  char * u = NULL;
  double scale = 1;
  float lc[3] = {0};
  float lw = 1.;
  int level = -1;
  Params params[] = {
    {"u", pstring, &u},
    {"scale", pdouble, &scale},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {"level", pint, &level},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!vectors (u,scale,lc,lw,level))
    return false;
}
/**
 * @brief Renders colored squares based on a scalar field or expression.
 *
 * Draws squares using the specified color or z-value expression, with options for value range, color mapping, face and line colors, orientation, and transparency.
 *
 * @param color Optional color specification or expression.
 * @param z Optional z-value expression for coloring.
 * @param min Minimum value for color mapping.
 * @param max Maximum value for color mapping.
 * @param spread Value spread for color scaling.
 * @param linear If true, applies linear color scaling.
 * @param fc Face color as an RGB array.
 * @param lc Line color as an RGB array.
 * @param n Normal vector for square orientation.
 * @param alpha Transparency value.
 * @return false if parameter parsing or rendering fails.
 */
else if (!strcmp (s, "squares")) {
  char * color = NULL;
  char * z = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  bool expr = false;
  coord n = {0,0,1};
  double alpha = 0;
  Params params[] = {
    {"color", pstring, &color},
    {"z", pstring, &z},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"n", pdouble, &n, 3},
    {"alpha", pdouble, &alpha},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!squares (color,z,min,max,spread,linear,map,fc,lc,expr,n,alpha))
    return false;
}
else if (!strcmp (s, "box")) {
  bool notics = false;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"notics", pbool, &notics},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!box (notics,lc,lw))
    return false;
}
else if (!strcmp (s, "isosurface")) {
  char * f = NULL;
  double v = 0.;
  char * color = NULL;
  double min = 0;
  double max = 0;
  double spread = 0;
  bool linear = false;
  Colormap map = jet;
  float fc[3] = {0};
  float lc[3] = {0};
  float lw = 1;
  bool expr = false;
  Params params[] = {
    {"f", pstring, &f},
    {"v", pdouble, &v},
    {"color", pstring, &color},
    {"min", pdouble, &min},
    {"max", pdouble, &max},
    {"spread", pdouble, &spread},
    {"linear", pbool, &linear},
    {"fc", pfloat, fc, 3},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!isosurface (f,v,color,min,max,spread,linear,map,fc,lc,lw,expr))
    return false;
}
else if (!strcmp (s, "travelling")) {
  double start = 0;
  double end = 0;
  float tx = 0;
  float ty = 0;
  float quat[4] = {0};
  float fov = 0;
  Params params[] = {
    {"start", pdouble, &start},
    {"end", pdouble, &end},
    {"tx", pfloat, &tx},
    {"ty", pfloat, &ty},
    {"quat", pfloat, quat, 4},
    {"fov", pfloat, &fov},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  travelling (start,end,tx,ty,quat,fov);
}
else if (!strcmp (s, "draw_string")) {
  char * str = NULL;
  int pos = 0;
  float size = 40;
  float lc[3] = {0};
  float lw = 1;
  Params params[] = {
    {"str", pstring, &str},
    {"pos", pint, &pos},
    {"size", pfloat, &size},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!draw_string (str,pos,size,lc,lw))
    return false;
}
/**
 * @brief Renders labels for a specified field with given line color and width.
 *
 * @param f Name of the field or label set to render.
 * @param lc RGB color array for the label outlines.
 * @param lw Line width for the label outlines.
 * @return false if parameter parsing or label rendering fails, true otherwise.
 */
else if (!strcmp (s, "labels")) {
  char * f = NULL;
  float lc[3] = {0};
  float lw = 1;
  Params params[] = {
    {"f", pstring, &f},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!labels (f,lc,lw))
    return false;
}
/**
 * @brief Draws lines from a file with specified color and width.
 *
 * Parses the file path, line color, and line width parameters, then renders the lines described in the file using the given visual properties.
 *
 * @return false if parameter parsing or line rendering fails.
 */
else if (!strcmp (s, "lines")) {
  char * file = NULL;
  float lc[3] = {0};
  float lw = 1.;
  Params params[] = {
    {"file", pstring, &file},
    {"lc", pfloat, lc, 3},
    {"lw", pfloat, &lw},
    {NULL}
  };
  if (!parse_params (params))
    return false;
  if (!lines (file,lc,lw))
    return false;
}
