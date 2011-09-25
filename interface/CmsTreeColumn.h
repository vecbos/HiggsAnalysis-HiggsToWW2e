// -*- C++ -*-
//---------------------------------------------------------------------------------
//
// $Id: CmsTreeColumn.h,v 1.2 2010/08/16 16:03:11 emanuele Exp $
//
// Description:
//    Class CmsTreeColumn
//    Nested class hierarchy to hold information about CmsTree columns.
// Author:
//    Emanuele Di Marco        (who stole a lot of code from BaBar RooTuple)
//
//---------------------------------------------------------------------------------

#ifndef CmsTreeColumn_h
#define CmsTreeColumn_h

#include <vector>
#include <string>

class TTree;
class TBranch;

using namespace std;

// Parent class (abstract):
class CmsTreeColumn {
public:
  CmsTreeColumn( const char* l ) : 
    label( l ), useDefValue( false ), pointer( 0 ), brp( 0 ) {}
  virtual ~CmsTreeColumn() {}
  virtual const std::string & getLabel() const { return label; }
  virtual TBranch* getBrPointer() { return brp; }
  virtual void* getPointer() { return pointer; }
  virtual void setPointer( void* p ) { pointer= p; }
  virtual void setUseDefValue( bool b ) { useDefValue= b; }
  virtual const bool & getUseDefValue() const { return useDefValue; }
  virtual void setDefValue() = 0;
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 ) = 0;
protected:
  std::string label; 
  bool useDefValue;
  void *pointer;
  TBranch* brp;
};

// Classes for bool:
class BoolCmsTreeColumn : public CmsTreeColumn {
public:
  BoolCmsTreeColumn( const char*, const bool &, const bool &, TTree* );
  virtual ~BoolCmsTreeColumn() { delete (char*)pointer; }
  virtual void setDefValue() { *(char*)pointer= defValue; }
  virtual void setValue( const void* p, CmsTreeColumn* cp= 0 ) { 
    *(char*)pointer= *(const bool*)p;
  }
private:
  bool defValue;
};
class BoolArrCmsTreeColumn : public CmsTreeColumn {
public:
  BoolArrCmsTreeColumn( const char*, const vector<bool> &, const bool &, 
		 TTree* );
  virtual ~BoolArrCmsTreeColumn() { delete[] (char*)pointer; }
  virtual void setDefValue() { 
    for( int i= 0; i < nmax; ++i ) ((char*)pointer)[i]= defValue; 
  }
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  bool defValue;
  int nmax;
};
class BoolDynArrCmsTreeColumn : public CmsTreeColumn {
public:
  BoolDynArrCmsTreeColumn( const char*, const vector<bool> &,
		    const bool &, CmsTreeColumn*, TTree* );
  virtual ~BoolDynArrCmsTreeColumn() { delete[] (char*)pointer; }
  virtual void setDefValue();
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  bool defValue;
  CmsTreeColumn* indexp;
};

// Classes for int:
class IntCmsTreeColumn : public CmsTreeColumn {
public:
  IntCmsTreeColumn( const char*, const int &, const int &, TTree* );
  virtual ~IntCmsTreeColumn() { delete (int*)pointer; }
  virtual void setDefValue() { *(int*)pointer= defValue; }
  virtual void setValue( const void* p, CmsTreeColumn* cp= 0 ) { 
    *(int*)pointer= *(const int*)p;
  }
private:
  int defValue;
};
class IntArrCmsTreeColumn : public CmsTreeColumn {
public:
  IntArrCmsTreeColumn( const char*, const vector<int> &, const int &, 
		TTree* );
  virtual ~IntArrCmsTreeColumn() { delete[] (int*)pointer; }
  virtual void setDefValue() { 
    for( int i= 0; i < nmax; ++i ) ((int*)pointer)[i]= defValue; 
  }
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  int defValue;
  int nmax;
};
class IntDynArrCmsTreeColumn : public CmsTreeColumn {
public:
  IntDynArrCmsTreeColumn( const char*, const vector<int> &,
                          const int &, CmsTreeColumn*, TTree* );
  virtual ~IntDynArrCmsTreeColumn() { delete[] (int*)pointer; }
  virtual void setDefValue();
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  int defValue;
  CmsTreeColumn* indexp;
};

// Classes for long:
class LongCmsTreeColumn : public CmsTreeColumn {
public:
  LongCmsTreeColumn( const char*, const uint64_t &, const uint64_t &, TTree* );
  virtual ~LongCmsTreeColumn() { delete (uint64_t*)pointer; }
  virtual void setDefValue() { *(uint64_t*)pointer= defValue; }
  virtual void setValue( const void* p, CmsTreeColumn* cp= 0 ) { 
    *(uint64_t*)pointer= *(const uint64_t*)p;
  }
private:
  uint64_t defValue;
};
class LongArrCmsTreeColumn : public CmsTreeColumn {
public:
  LongArrCmsTreeColumn( const char*, const vector<uint64_t> &, const uint64_t &, 
                        TTree* );
  virtual ~LongArrCmsTreeColumn() { delete[] (uint64_t*)pointer; }
  virtual void setDefValue() { 
    for( int i= 0; i < nmax; ++i ) ((uint64_t*)pointer)[i]= defValue; 
  }
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  uint64_t defValue;
  int nmax;
};
class LongDynArrCmsTreeColumn : public CmsTreeColumn {
public:
  LongDynArrCmsTreeColumn( const char*, const vector<uint64_t> &,
                           const uint64_t &, CmsTreeColumn*, TTree* );
  virtual ~LongDynArrCmsTreeColumn() { delete[] (uint64_t*)pointer; }
  virtual void setDefValue();
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  uint64_t defValue;
  CmsTreeColumn* indexp;
};

// Classes for float:
class FloatCmsTreeColumn : public CmsTreeColumn {
public:
  FloatCmsTreeColumn( const char*, const float &, const float &, TTree* );
  virtual ~FloatCmsTreeColumn() { delete (float*)pointer; }
  virtual void setDefValue() { *(float*)pointer= defValue; }
  virtual void setValue( const void* p, CmsTreeColumn* cp= 0 ) { 
    *(float*)pointer= *(const float*)p; 
  }
private:
  float defValue;
};
class FloatArrCmsTreeColumn : public CmsTreeColumn {
public:
  FloatArrCmsTreeColumn( const char*, const vector<float> &, const float &, 
		  TTree* );
  virtual ~FloatArrCmsTreeColumn() { delete[] (float*)pointer; }
  virtual void setDefValue() { 
    for( int i= 0; i < nmax; ++i ) ((float*)pointer)[i]= defValue; 
  }
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  float defValue;
  int nmax;
};
class FloatDynArrCmsTreeColumn : public CmsTreeColumn {
public:
  FloatDynArrCmsTreeColumn( const char*, const vector<float> &,
		     const float &, CmsTreeColumn*, TTree* );
  virtual ~FloatDynArrCmsTreeColumn() { delete[] (float*)pointer; }
  virtual void setDefValue();
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  float defValue;
  CmsTreeColumn* indexp;
};

// Classes for double:
class DoubleCmsTreeColumn : public CmsTreeColumn {
public:
  DoubleCmsTreeColumn( const char*, const double &, const double &, TTree* );
  virtual ~DoubleCmsTreeColumn() { delete (double*)pointer; }
  virtual void setDefValue() { *(double*)pointer= defValue; }
  virtual void setValue( const void* p, CmsTreeColumn* cp= 0 ) { 
    *(double*)pointer= *(const double*)p; 
  }
private:
  double defValue;
};
class DoubleArrCmsTreeColumn : public CmsTreeColumn {
public:
  DoubleArrCmsTreeColumn( const char*, const vector<double> &, 
		   const double &, TTree* );
  virtual ~DoubleArrCmsTreeColumn() { delete[] (double*)pointer; }
  virtual void setDefValue() { 
    for( int i= 0; i < nmax; ++i ) ((double*)pointer)[i]= defValue; 
  }
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  double defValue;
  int nmax;
};
class DoubleDynArrCmsTreeColumn : public CmsTreeColumn {
public:
  DoubleDynArrCmsTreeColumn( const char*, const vector<double> &,
		      const double &, CmsTreeColumn*, TTree* );
  virtual ~DoubleDynArrCmsTreeColumn() { delete[] (double*)pointer; }
  virtual void setDefValue();
  virtual void setValue( const void*, CmsTreeColumn* cp= 0 );
private:
  double defValue;
  CmsTreeColumn* indexp;
};

// Classes for string:
class StringCmsTreeColumn : public CmsTreeColumn {
public:
  StringCmsTreeColumn( const char*, const string &, const string &, TTree* );
  virtual ~StringCmsTreeColumn() { delete[] (char*)pointer; }
  virtual void setDefValue();
  virtual void setValue( const void* p, CmsTreeColumn* cp= 0 );
private:
  string defValue;
};

#endif  // CmsTreeColumn_h
