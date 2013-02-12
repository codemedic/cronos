#ifndef RECT_HPP
#define RECT_HPP

class TRect {
public:
  int left,top,right,bottom;
  TRect(int ix,int iy,int iw,int ih) : left(ix),top(iy),right(ix+iw-1),bottom(iy+ih-1) {};
  TRect() {};

  int Height() {return (bottom-top+1);}
  int Width() {return (right-left+1);}
  TRect& operator=(const TRect& other)
  {
    left = other.left;   right = other.right;
    top = other.top;  bottom = other.bottom;
    return *this;
  }
};

#endif

