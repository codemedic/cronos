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
 * information about a time series model.
 *
 */

#include "includes.h"
#include "tsmrender.h"

TSMRenderer::TSMRenderer(TimeSeriesModel *tsm, TimeSeries *ts,
			 vector<Transformation*> *tl)
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

  greektag = Gtk::TextBuffer::Tag::create("greektag");
  //greektag->property_family() = "symbol";
  tagtable->add(greektag);

  Update(tsm,ts,tl);
}

void TSMRenderer::on_realize()
{
  TextView::on_realize();
}

void TSMRenderer::Update(TimeSeriesModel* tsm, TimeSeries* ts,
			 vector<Transformation *> *transformlist)
{
  Transformation *trans;
  ostringstream os;

  if (tsm!=NULL)
    {
      tsm->RenderInto(os);
      os << endl << endl;
    }

  if (ts!=NULL)
    {
      ts->RenderInto(os);
      os << endl;
    }

  // now render the transformations
  if (transformlist!=NULL)
    if (!transformlist->empty()) {
      os << "Transforms:" << endl;
      int i,n=transformlist->size();
      for (i=0 ; i<n ; ++i)
	{
	  os << "   ";
	  (*transformlist)[i]->DescribeIn(os);
	  os << endl;
	}
      os << endl;
    }

  os << "\\n";
  os << ends;


  // now go through and insert the string
  string bigstring = os.str();
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
	      currenttag = &supertag; // super/sub? script?
	      break;
	    case 'g':
	      currenttag = &greektag;
	      break;
	    }
	  // delete first two chars
	  bigstring.erase(0,2);
	}
      // find next occurrence of backslash
      pos = bigstring.find('\\');
      if (pos==-1)
	pos = bigstring.size();

 
      if (currenttag!=NULL)
	tbuf->insert_with_tag(tbuf->end(), bigstring.substr(0,pos), *currenttag);
      else
	tbuf->insert(tbuf->end(), bigstring.substr(0,pos));
      bigstring.erase(0,pos);
    }

  set_buffer(tbuf);
}
