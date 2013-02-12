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
 * This module defines a class which displays
 * information and can print bold, greek, superscript, etc.
 *
 */

#include "includes.h"
#include "textrender.h"

TextRenderer::TextRenderer()
  : Gtk::TextView()
{
  set_editable(false);
  set_wrap_mode(Gtk::WRAP_WORD);
  set_name("MITSMTextWidget");

  tbuf = Gtk::TextBuffer::create();
  tagtable = tbuf->get_tag_table();
  boldtag = Gtk::TextBuffer::Tag::create("bolding");
  boldtag->property_weight() = Pango::WEIGHT_BOLD;
  boldtag->property_scale() = 1.2;
  tagtable->add(boldtag);

  supertag = Gtk::TextBuffer::Tag::create("supertag");
  supertag->property_rise() = 2500;
  supertag->property_scale() = 0.7;
  tagtable->add(supertag);

  subtag = Gtk::TextBuffer::Tag::create("subtag");
  subtag->property_rise() = -2500;
  subtag->property_scale() = 0.7;
  tagtable->add(subtag);

  tagmono = Gtk::TextBuffer::Tag::create("tagmono");
  tagmono->property_family()="monospace";
  tagtable->add(tagmono);

  tag0 = Gtk::TextBuffer::Tag::create("tag0");
  tag0->property_foreground() = "#e00000";
  tag0->property_family()="monospace";
  tagtable->add(tag0);

  tag1 = Gtk::TextBuffer::Tag::create("tag1");
  tag1->property_background() = "#e06060";
  tag1->property_family()="monospace";
  tagtable->add(tag1);

  tag2 = Gtk::TextBuffer::Tag::create("tag2");
  tag2->property_background() = "#e0c0c0";
  tag2->property_family()="monospace";
  tagtable->add(tag2);

  tag3 = Gtk::TextBuffer::Tag::create("tag3");
  tag3->property_background() = "#c0c0c0";
  tag3->property_family()="monospace";
  tagtable->add(tag3);

  tag4 = Gtk::TextBuffer::Tag::create("tag4");
  tag4->property_background() = "#60c060";
  tag4->property_family()="monospace";
  tagtable->add(tag4);

  tag5 = Gtk::TextBuffer::Tag::create("tag5");
  tag5->property_family()="monospace";
  tag5->property_background() = "#00e000";
  tagtable->add(tag5);

  tagunderline = Gtk::TextBuffer::Tag::create("tagunderline");
  tagunderline->property_underline() = Pango::UNDERLINE_SINGLE;
  tagtable->add(tagunderline);
}

void TextRenderer::on_realize()
{
  TextView::on_realize();
}

void TextRenderer::Update(const char *temps)
{
  // now go through and insert the string
  string bigstring(temps);
  int pos;
  Glib::RefPtr<Gtk::TextBuffer::Tag> *currenttag = NULL;
  tbuf->erase(tbuf->begin(), tbuf->end());

  while (bigstring.size()>0)
    {
      // pick out first bit
      if (bigstring[0]=='\\')
	{
	  switch (bigstring[1])
	    {
	    case 'n':
	      currenttag = NULL;
	      break;
	    case 'b': // bold for line
	      currenttag = &boldtag;
	      break;
	    case 's':
	      currenttag = &supertag; // superscript
	      break;
	    case 'd':
	      currenttag = &subtag;   // subscript
	      break;
	    case '0':
	      currenttag = &tag0;
	      break;
	    case '1':
	      currenttag = &tag1;
	      break;
	    case '2':
	      currenttag = &tag2;
	      break;
	    case '3':
	      currenttag = &tag3;
	      break;
	    case '4':
	      currenttag = &tag4;
	      break;
	    case '5':
	      currenttag = &tag5;
	      break;
	    case 'm':
	      currenttag = &tagmono;
	      break;
	    case 'u':
	      currenttag = &tagunderline;
	      break;
	    }
	  // delete first two chars
	  bigstring.erase(0,2);
	}
      // find next occurrence of backslash
      pos = bigstring.find('\\');
      if (pos==-1)
	pos = bigstring.size();


      try 
	{
	  if (currenttag!=NULL)
	    tbuf->insert_with_tag(tbuf->end(), bigstring.substr(0,pos), *currenttag);
	  else
	    tbuf->insert(tbuf->end(), bigstring.substr(0,pos));
	}
      catch(char *str)
	{
	  cout << "Exception:" << str << endl;
	}
      bigstring.erase(0,pos);
    }

  set_buffer(tbuf);
}
