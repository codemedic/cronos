/*
 * -------------------------------------------------------------------
 *
 * Copyright 2004 Anthony Brockwell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * -------------------------------------------------------------------
 *
 * This module contains a class which writes to encapsulated
 * postscript files, treating them like standard graphics
 * devices, with calls to generate lines, rectangles, etc.
 * March, 2001, Anthony Brockwell.
 *
 */


#include "epsdc.h"
#include "rect.h"
#include <iostream.h>
#include <fstream.h>

TEpsPlot::TEpsPlot(TRect boundingbox, const char *filename)
{
  outfile = new ofstream(filename);
  topy = boundingbox.Height();
  bottomy = 0;
  scale = 576.0/boundingbox.Width();  // 8 inches wide

  // write header
  *outfile << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
  *outfile << "%%Creator: Brockwell Software" << endl;
  *outfile << "%%Title: test1.eps" << endl;
  *outfile << "%%DocumentNeededFonts: Times-Roman" << endl;
  *outfile << "%%LanguageLevel: 2" << endl;
  *outfile << "%%Pages: 0" << endl;
  *outfile << "%%Orientation: Portrait" << endl;
  *outfile << "%%BoundingBox: 0 0 " << boundingbox.Width()*scale << " "
  	 << boundingbox.Height()*scale << endl;
  *outfile << "%%EndComments" << endl;
  *outfile << endl << endl;

  *outfile << "/Times-Roman findfont" << endl;
  *outfile << "30 scalefont" << endl;
  *outfile << "setfont" << endl;

  *outfile << endl;

  *outfile << "/L {lineto} bind def" << endl;
  *outfile << "/M {moveto} bind def" << endl;
  *outfile << "/S {stroke} bind def" << endl;
  *outfile << "/N {newpath} bind def" << endl;
  *outfile << "/D {setdash} bind def" << endl;

  *outfile << scale << " " << scale << " scale" << endl;
  *outfile << "2 setlinewidth" << endl;

  inpath = false;
}

TEpsPlot::~TEpsPlot()
{
  FinishPath();
  *outfile << endl << "showpage" << endl << "%%Trailer" << endl
    << "%%EOF" << endl;
  outfile->close();
}

int TEpsPlot::tof(int y)
{
  return topy - y;
}

void TEpsPlot::TextOut(int x, int y, const char *cp, bool center)
{
  FinishPath();
  *outfile << x << " (" << cp << ") stringwidth pop ";
  if (center)
    *outfile << "2 div ";
  *outfile << "sub ";
  *outfile << tof(y);
  *outfile <<  " M" << endl;
  *outfile << "(" << cp << ") show" << endl;
}

void TEpsPlot::FinishPath()
{
  if (inpath) {
    *outfile << "S" << endl;  // finish previous path
    inpath = false;
  }
}

void TEpsPlot::SetLineMode(int mode, int width)
{
  FinishPath();
  switch (mode)
    {
      case 1:
        *outfile << "[4 10] 0 D" << endl;
        break;
      default:
        *outfile << "[] 0 D" << endl;
        break;
    }
  *outfile << width << " setlinewidth" << endl;
}

void TEpsPlot::FillRect(int x1, int y1, int x2, int y2, 
  int r1, int g1, int b1)
{
  float r=r1/256.0,g=g1/256.0,b=b1/256.0;
  DrawRect(x1,y1,x2,y2);
  *outfile << r << " " << g << " " << b << " setrgbcolor" << endl;
  *outfile << " fill" << endl;
  *outfile << "0 0 0 setrgbcolor" << endl;
}

void TEpsPlot::DrawRect(int x1, int y1, int x2, int y2)
{
  MoveTo(x1,y1);
  LineTo(x2,y1);
  LineTo(x2,y2);
  LineTo(x1,y2);
  LineTo(x1,y1);
}

void TEpsPlot::SetFontSize(int w, int h)
{
  FinishPath();
  *outfile << "/Times-Roman findfont" << endl
    << "[" << h << " 0 0 " << h << " 0 0] makefont" << endl
  	 << "setfont" << endl;
}

bool TEpsPlot::MoveTo(int x, int y)
{
  FinishPath();
  *outfile << "N" << endl << x << " " << tof(y) << " M" << endl;
  return true;
}

bool TEpsPlot::LineTo(int x, int y)
{
  *outfile << x << " " << tof(y) << " L" << endl;
  inpath = true;
  return true;
}



