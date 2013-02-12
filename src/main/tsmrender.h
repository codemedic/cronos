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

#ifndef TSMRENDER_HPP
#define TSMRENDER_HPP

class TSMRenderer : public Gtk::TextView {
protected:
  TimeSeriesModel *tsmp;

  Glib::RefPtr<Gtk::TextBuffer::TagTable> tagtable;
  Glib::RefPtr<Gtk::TextBuffer::Tag> boldtag, supertag, greektag;
  Glib::RefPtr<Gtk::TextBuffer> tbuf;

  /*
  Gdk::Color white,black,grey,red;
  Gdk::GC gc;
  Gdk::Window window;
  Gdk::Font *f1, *fb, *fs;
  */

  void on_realize();

public:
  TSMRenderer(TimeSeriesModel *, TimeSeries *, vector<Transformation*> *tl);
  void Update(TimeSeriesModel *, TimeSeries *, vector<Transformation*> *tl);
};

#endif

