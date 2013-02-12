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
 */

#ifndef EPSDC_HPP
#define EPSDC_HPP

#include "rect.h"
#include <fstream.h>

class TEpsPlot
{
protected:
  ofstream *outfile;
  bool inpath;
  int topy,bottomy;
  double scale;
  int tof(int);
  
public:
  TEpsPlot(TRect bbox, const char *filename);
  ~TEpsPlot();
  
  // things I have taken over
  void FinishPath();
  bool MoveTo(int,int);
  bool LineTo(int,int);
  void FillRect(int x1, int y1, int x2, int y2, int r,int g,int b);
  void FillRect(TRect r1, int r, int g, int b)
  { FillRect(r1.left, r1.top, r1.right, r1.bottom, r, g, b); };
  void DrawRect(int x1, int y1, int x2, int y2);
  void DrawRect(TRect r1)
  { DrawRect(r1.left, r1.top, r1.right, r1.bottom);};
  void TextOut(int,int,const char *,bool center);
  
  void SetLineMode(int md, int w=2);
  void SetFontSize(int w, int h);
};

#endif
