// /*
//   This file provides some tools for drawing to post-script files
// */
// 
// #ifndef PS_H__
// #define PS_H__
// 
// /*
//   This just opens and writes a header to our PS file
// */
// FILE* PS_Start(char *fn, float w, float h)
// {
//   FILE *f;
// 
//   // Open file
//   if( !(f = fopen(fn,"w")) )
//     return NULL;
// 
//   // Write header (procedure definitions)
//   fprintf(f,
// 	  "%%!PS-Adobe-3.0 EPSF-3.0\n"
// 	  "%%%%Creator: Created by ps.h (Aaron Vose)\n"
// 	  "%%%%BoundingBox: 0 0 %d %d\n"
// 	  "%%%%LanguageLevel: 2\n"
// 	  "%%%%Pages: 1\n"
// 	  "%%%%DocumentData: Clean7Bit\n\n"
// 	  "/inch {72 mul} def\n"
// 	  "/flip {%f inch sub -1 mul} def\n\n"
// 	  "save\n"
// 	  "%f inch %f inch translate\n"
// 	  "180 rotate\n\n",
// 	  (int)(w*72),(int)(h*72),
// 	  w,
// 	  w,h
// 	  );
// 
//   // Just return the file pointer
//   return f;
// }
// 
// /*
//   Changes the current color
// */
// void PS_Color(FILE *f, float r, float g, float b)
// {
//   // Set the current color
//   fprintf(f,"%f %f %f setrgbcolor\n",r,g,b);
// }
// 
// /*
//   Changes the current color (grayscale)
// */
// void PS_Gray(FILE *f, float i)
// {
//   // Set the current color (grayscale)
//   fprintf(f,"%f setgray\n",i);
// }
// 
// /*
//   Rotates the 'page' a degrees counter-clockwise
// */
// void PS_Rotate(FILE *f, float a)
// {
//   // Set the current color (grayscale)
//   fprintf(f,"%f rotate\n",a);
// }
// 
// /*
//   Draws a line
// */
// void PS_Line(FILE *f, float x1, float y1, float x2, float y2)
// {
//   // Define the line path
//   fprintf(f,
// 	  "newpath\n"
// 	  ".05 setlinewidth\n"
// 	  "%f inch flip %f inch moveto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "closepath\n",
// 	  x1,y1,
// 	  x2,y2
// 	  );
// 
//   // Stroke or fill the path
//   fprintf(f, "stroke\n");
// }
// 
// /*
//   Draws a quad
// */
// void PS_Quad(FILE *f, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4, int fill)
// {
//   // Define the box path
//   fprintf(f,
// 	  "newpath\n"
// 	  "%f inch flip %f inch moveto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "closepath\n",
// 	  x1,y1,
// 	  x2,y2,
// 	  x3,y3,
// 	  x4,y4
// 	  );
// 
//   // Stroke or fill the path
//   if(fill) fprintf(f, "fill\n");
//   else     fprintf(f, "stroke\n");
// }
// 
// /*
//   Draws a box
// */
// void PS_Box(FILE *f, float x, float y, float w, float h, int fill)
// {
//   // Define the box path
//   fprintf(f,
// 	  "newpath\n"
// 	  "%f inch flip %f inch moveto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "%f inch flip %f inch lineto\n"
// 	  "closepath\n",
// 	  x,y,
// 	  x+w,y,
// 	  x+w,y+h,
// 	  x,y+h
// 	  );
// 
//   // Stroke or fill the path
//   if(fill) fprintf(f, "fill\n");
//   else     fprintf(f, "stroke\n");
// }
// 
// /*
//   Draws a filled or outlined circle
// */
// void PS_Circle(FILE *f, float x, float y, float r, int fill)
// {
//   // Draw a circle to the path
//   fprintf(f,
// 	  "newpath\n"
// 	  "%f inch flip %f inch %f inch 0 360 arc\n"
// 	  "closepath\n",
// 	  x,y,r);
//   
//   // Stroke or fill the path
//   if(fill) fprintf(f, "fill\n");
//   else     fprintf(f, "stroke\n");
// }
// 
// /*
//   Draws text at x,y with size x (times font)
// */
// void PS_Text(FILE *f, float x, float y, int s, char *t)
// {
// 
//   // Set up the font
//   fprintf(f,
// 	  "/Times-Roman findfont\n"
// 	  "%d scalefont\n"
// 	  "setfont\n",
// 	  s
// 	  );
// 
//   // Draw the string
//   fprintf(f,
// 	  "newpath\n"
// 	  "%f inch flip %f inch moveto\n"
// 	  "(%s) show\n"
// 	  "closepath\n",
// 	  x,y,
// 	  t
// 	  );
// }
// 
// /*
//   Writes a few finishing things to the file and then closes it
// */
// void PS_End(FILE *f)
// {
//   // Write trailer and exit
//   fprintf(f,"\nrestore\n%%EOF\n");
//   fclose(f);
// }
// 
// 
// #endif

/*
  This file provides some tools for drawing to post-script files
*/

#ifndef PS_H__
#define PS_H__

/*
  This just opens and writes a header to our PS file
*/
FILE* PS_Start(char *fn, float w, float h)
{
  FILE *f;

  // Open file
  if( !(f = fopen(fn,"w")) )
    return NULL;

  // Write header (procedure definitions)
  fprintf(f,
	  "%%!PS-Adobe-3.0 EPSF-3.0\n"
	  "%%%%Creator: Created by ps.h (Aaron Vose)\n"
	  "%%%%BoundingBox: 0 0 %d %d\n"
	  "%%%%LanguageLevel: 2\n"
	  "%%%%Pages: 1\n"
	  "%%%%DocumentData: Clean7Bit\n\n"
	  "/inch {72 mul} def\n"
	  "/flip {%f inch sub -1 mul} def\n\n"
	  "save\n",
	  (int)(w*72),(int)(h*72),
	  h
	  );

  // Just return the file pointer
  return f;
}

/*
  Changes the current color
*/
void PS_Color(FILE *f, float r, float g, float b)
{
  // Set the current color
  fprintf(f,"%f %f %f setrgbcolor\n",r,g,b);
}

/*
  Changes the current color (grayscale)
*/
void PS_Gray(FILE *f, float i)
{
  // Set the current color (grayscale)
  fprintf(f,"%f setgray\n",i);
}

/*
  Rotates the 'page' a degrees counter-clockwise
*/
void PS_Rotate(FILE *f, float a)
{
  // Set the current color (grayscale)
  fprintf(f,"%f rotate\n",a);
}

void PS_Translate(FILE *f, float x, float y)
{
  fprintf(f,"%f inch %f inch translate\n",x,y);
}

/*
  Draws a line
*/
void PS_Line(FILE *f, float x1, float y1, float x2, float y2)
{
  // Define the line path
  fprintf(f,
	  "newpath\n"
	  ".05 setlinewidth\n"
	  "%f inch %f inch flip moveto\n"
	  "%f inch %f inch flip lineto\n"
	  "closepath\n",
	  x1,y1,
	  x2,y2
	  );

  // Stroke or fill the path
  fprintf(f, "stroke\n");
}

/*
  Draws a quad
*/
void PS_Quad(FILE *f, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4, int fill)
{
  // Define the box path
  fprintf(f,
	  "newpath\n"
	  "%f inch %f inch flip moveto\n"
	  "%f inch %f inch flip lineto\n"
	  "%f inch %f inch flip lineto\n"
	  "%f inch %f inch flip lineto\n"
	  "closepath\n",
	  x1,y1,
	  x2,y2,
	  x3,y3,
	  x4,y4
	  );

  // Stroke or fill the path
  if(fill) fprintf(f, "fill\n");
  else     fprintf(f, "stroke\n");
}

/*
  Draws a box
*/
void PS_Box(FILE *f, float x, float y, float w, float h, int fill)
{
  // Define the box path
  fprintf(f,
	  "newpath\n"
	  "%f inch %f inch flip moveto\n"
	  "%f inch %f inch flip lineto\n"
	  "%f inch %f inch flip lineto\n"
	  "%f inch %f inch flip lineto\n"
	  "closepath\n",
	  x,y,
	  x+w,y,
	  x+w,y+h,
	  x,y+h
	  );

  // Stroke or fill the path
  if(fill) fprintf(f, "fill\n");
  else     fprintf(f, "stroke\n");
}

/*
  Draws a filled or outlined circle
*/
void PS_Circle(FILE *f, float x, float y, float r, int fill)
{
  // Draw a circle to the path
  fprintf(f,
	  "newpath\n"
	  "%f inch %f inch flip %f inch 0 360 arc\n"
	  "closepath\n",
	  x,y,r);
  
  // Stroke or fill the path
  if(fill) fprintf(f, "fill\n");
  else     fprintf(f, "stroke\n");
}

/*
  Draws text at x,y with size x (times font)
*/
void PS_Text(FILE *f, float x, float y, int s, char *t)
{

  // Set up the font
  fprintf(f,
	  "/Times-Roman findfont\n"
	  "%d scalefont\n"
	  "setfont\n",
	  s
	  );

  // Draw the string
  fprintf(f,
	  "newpath\n"
	  "%f inch %f inch flip moveto\n"
	  "(%s) show\n"
	  "closepath\n",
	  x,y,
	  t
	  );
}

/*
  Writes a few finishing things to the file and then closes it
*/
void PS_End(FILE *f)
{
  // Write trailer and exit
  fprintf(f,"\nrestore\n%%EOF\n");
  fclose(f);
}


#endif
