// -*- C++ -*-
//-----------------------------------------------------------------------
//
// Description:
//    Class CmsTree
//    see the CmsTree.h for documentation
// Author:
//    Emanuele Di Marco        (who stole a lot of code from BaBar RooTuple)
//
//---------------------------------------------------------------------------------


#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTree.h"

#include <iostream>
#include <string.h>

#include "TTree.h"
#include "TBranch.h"

using namespace std;


// Constructor to create a tuple with name and title:
CmsTree::CmsTree( const char* name, const char* title ) {

  treep = new TTree( name, title );
  
}

// Destructor:
CmsTree::~CmsTree() {

  delete treep;

}

// Column booking/filling. All these have the same name - column(...)

// Specify the data for a column. The string is to the key to
// the column, so it must be unique.


// ====== Bool type ======
// Make/fill column with a single value
void CmsTree::column( const char* label, 
		       bool value,
		       bool defval,
		       const char* block ) {
  
  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &value );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new column:
    colp= new BoolCmsTreeColumn( label, value, defval, treep );
    lpmap.append( colp );
  }

}
// Make/fill column-array. Length is fixed at creation time.
void CmsTree::column( const char* label, 
		       const vector<bool>& v, 
		       bool defval,
		       const char* block ){

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new Column:
    colp= new BoolArrCmsTreeColumn( label, v, defval, treep );
    lpmap.append( colp );
  }

}
// Make/fill column-array. Length is variable and is taken from 
// another column.
void CmsTree::column( const char* label, 
		       const vector<bool>& v, 
		       const char* ilab,
		       bool defval,
		       const char* block ) {
  
  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v, lpmap.getColumn( ilab ) );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new branch:
    CmsTreeColumn* indexp= lpmap.getColumn( ilab );
    colp= new BoolDynArrCmsTreeColumn( label, v, defval, indexp, treep );
    lpmap.append( colp );
  }
  
}




// ====== Int type ======
// Make/fill column with a single value
void CmsTree::column( const char* label, 
		       int value, 
		       int defval, 
		       const char* block) {
  
  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &value );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new column:
    colp= new IntCmsTreeColumn( label, value, defval, treep );
    lpmap.append( colp );
  }

}
// Make/fill column-array. Length is fixed at creation time.
void CmsTree::column( const char* label, 
		       const vector<int> & v, 
		       int defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new Column:
    colp= new IntArrCmsTreeColumn( label, v, defval, treep );
    lpmap.append( colp );
  }

}
// Make/fill column-array. Length is variable and is taken from 
// another column.
void CmsTree::column( const char* label, 
		       const vector<int> & v, 
		       const char* ilab,
		       int defval,
		       const char* block) {
  
  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v, lpmap.getColumn( ilab ) );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new branch:
    CmsTreeColumn* indexp= lpmap.getColumn( ilab );
    colp= new IntDynArrCmsTreeColumn( label, v, defval, indexp, treep );
    lpmap.append( colp );
  }
  
}















// ====== Float type ======
// Make/fill column with a single value
void CmsTree::column( const char* label, 
		       float value,
		       float defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &value );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new column:
    colp= new FloatCmsTreeColumn( label, value, defval, treep );
    lpmap.append( colp );
  }


}
// Make/fill column-array. Length is fixed at creation time.
void CmsTree::column( const char* label, 
		       const vector<float> & v, 
		       float defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new Column:
    colp= new FloatArrCmsTreeColumn( label, v, defval, treep );
    lpmap.append( colp );
  }


}
// Make/fill column-array. Length is variable and is taken from 
// another column.
void CmsTree::column( const char* label, 
		       const vector<float> & v, 
		       const char *ilab,
		       float defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v, lpmap.getColumn( ilab ) );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new branch:
    CmsTreeColumn* indexp= lpmap.getColumn( ilab );
    colp= new FloatDynArrCmsTreeColumn( label, v, defval, indexp, treep );
    lpmap.append( colp );
  }

}

// ====== Double type ======
// Make/fill column with a single value
void CmsTree::column( const char* label, 
		       double value,
		       double defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &value );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new column:
    colp= new DoubleCmsTreeColumn( label, value, defval, treep );
    lpmap.append( colp );
  }


}
// Make/fill column-array. Length is fixed at creation time.
void CmsTree::column( const char* label, 
		       const vector<double> & v, 
		       double defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new Column:
    colp= new DoubleArrCmsTreeColumn( label, v, defval, treep );
    lpmap.append( colp );
  }


}
// Make/fill column-array. Length is variable and is taken from 
// another column.
void CmsTree::column( const char* label, 
		       const vector<double> & v, 
		       const char *ilab,
		       double defval,
		       const char* block) {

  CmsTreeColumn* colp= lpmap.getColumn( label );
  if( colp ) {
    // Column exists, fill corresponding memory location with value:
    colp->setValue( &v, lpmap.getColumn( ilab ) );
    colp->setUseDefValue( false );
  }
  else {
    // Create a new branch:
    CmsTreeColumn* indexp= lpmap.getColumn( ilab );
    colp= new DoubleDynArrCmsTreeColumn( label, v, defval, indexp, treep );
    lpmap.append( colp );
  }

}



// Dump all the data into the ntuple and then clear
int CmsTree::dumpData() { 

  vector<CmsTreeColumn*> clist= lpmap.getMap();
  for( int i= 0; i < (int)clist.size(); ++i ) {
    if( clist[i]->getUseDefValue() ) clist[i]->setDefValue();
    clist[i]->setUseDefValue( true );
  }
  int row_number = (int) treep->GetEntries();
  treep->Fill();
  return row_number;
}

// Set all the data to their default values:
void CmsTree::clearData() { 

  vector<CmsTreeColumn*> clist= lpmap.getMap();
  for( int i= 0; i < (int)clist.size(); ++i ) {
    clist[i]->setDefValue();
    clist[i]->setUseDefValue( true );
  }
  return; 

}

// Return title of ntuple:
std::string CmsTree::title() const { return treep->GetTitle(); }

// number of columns
int CmsTree::nColumns() const {

  vector<CmsTreeColumn*> map= lpmap.getMap();
  return map.size();

}

// Return label for a particular column with index i:
std::string CmsTree::label( int i ) const {

  vector<CmsTreeColumn*> map= lpmap.getMap();
  string str;
  if( i >= 0 && i < (int)map.size() ) str= map[i]->getLabel();
  else str= "unknown column index";
  return str;

}

// Print info about ntuple:
void CmsTree::print( ostream & o ) const { 

  cout << "CmsTree: ntuple " << title() << " has " << nColumns() 
       << " columns." << endl;
  cout << "Complete printout follows: " << endl;
  treep->Print(); 


}

// Reset:
void CmsTree::reset() {

  treep->Reset();
  return;

}


// Member functions of nested class LPmap:

// Return pointer to column with label l, cache last returned column:
CmsTreeColumn* CmsTree::LPmap::getColumn( const char* l ) const {

  static string oldstr;
  static int i;
  if( oldstr.c_str() != l ) {
    oldstr= l;
    for( i= 0; i < (int)map.size(); ++i ) if( map[i]->getLabel() == l ) break;
  }
  CmsTreeColumn* colp= 0;
  if( i < (int)map.size() ) colp=map[i];
  return colp;

}

