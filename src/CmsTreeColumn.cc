// -*- C++ -*-
//---------------------------------------------------------------------------------
//
// $Id: CmsTreeColumn.cc,v 1.1 2007/04/13 13:05:16 emanuele Exp $
//
// Description:
//    Class CmsTreeColumn
//    Nested class hierarchy to hold information about CmsTree columns.
// Author:
//    Emanuele Di Marco        (who stole a lot of code from BaBar RooTuple)
//
//---------------------------------------------------------------------------------

#include <memory>

#include "HiggsAnalysis/HiggsToWW2e/interface/CmsTreeColumn.h"
#include <TTree.h>
#include <TBranch.h>
#include <iostream>

using namespace std;


// Bool columns, stored as char (= 8 bit signed int):
BoolCmsTreeColumn::BoolCmsTreeColumn( const char* l, const bool & v, const bool & d, 
			TTree* tp ) : 
  CmsTreeColumn( l ), defValue( d ) {
				 
  // Create a new branch:
  pointer= new char;
  *(char*)pointer= v;
  std::string leafs( l ) ;
  leafs+= "/b";
  brp= tp->Branch( label.c_str(), pointer, leafs.c_str(), 8000 );

}
BoolArrCmsTreeColumn::BoolArrCmsTreeColumn( const char* l, 
			      const vector<bool> & v, 
			      const bool & d, TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ) {

  // Create a new branch:
  nmax= v.size();
  char* bp= new char[nmax];
  pointer= bp;
  for( int i= 0; i < nmax; ++i ) bp[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  char buf[33];
  sprintf( buf, "%i", nmax );
  leafs+= buf;
  leafs+= "]/b";
  brp= tp->Branch( label.c_str(), &bp[0], leafs.c_str(), 8000 );

}
void BoolArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) { 

  const vector<bool>* vp= (const vector<bool>*) p;
  if( (int)vp->size() < nmax ) {
    std::cerr << "BoolArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    for( int i= 0; i < nmax; ++i ) ((char*)pointer)[i]= (*vp)[i];
  }

}
BoolDynArrCmsTreeColumn::BoolDynArrCmsTreeColumn( const char* l, 
				    const vector<bool> & v, 
				    const bool & d, CmsTreeColumn* inp,
				    TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ), indexp( inp ) {

  // Make a new branch:
  int* np= (int*) indexp->getPointer();
  char* bp= new char[*np];
  pointer= bp;
  for( int i= 0; i < *np; ++i ) bp[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  leafs+= indexp->getLabel();
  leafs+= "]/b";
  brp= tp->Branch( label.c_str(), &bp[0], leafs.c_str(), 8000 );

}
void BoolDynArrCmsTreeColumn::setDefValue() {

  if( pointer ) delete [] (char*)pointer;
  int nmax= *((int*)(indexp->getPointer()));
  char* bp= new char[nmax];
  pointer= bp;
  for( int i= 0; i < nmax; ++i ) bp[i]= defValue; 
  brp->SetAddress( &bp[0] );

}
void BoolDynArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) {

  const vector<bool>* vp= (const vector<bool>*) p;
  int* np= (int*) cp->getPointer();
  if( *np > (int)vp->size() ) {
    std::cerr << "BoolDynArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    char* bp= new char[*np];
    if( pointer ) delete[] (char*)pointer;
    pointer= bp;
    for( int i= 0; i < *np; ++i ) bp[i]= (*vp)[i];
    brp->SetAddress(  &bp[0] );
  }

}


// Int columns:
IntCmsTreeColumn::IntCmsTreeColumn( const char* l, const int & v, const int & d, 
		      TTree* tp ) : 
  CmsTreeColumn( l ), defValue( d ) {
				 
  // Create a new branch:
  pointer= new int;
  *(int*)pointer= v;
  std::string leafs( l ) ;
  leafs+= "/I";
  brp= tp->Branch( label.c_str(), pointer, leafs.c_str(), 8000 );

}
IntArrCmsTreeColumn::IntArrCmsTreeColumn( const char* l, 
			    const vector<int> & v, 
			    const int & d, TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ) {

  // Create a new branch:
  nmax= v.size();
  int* ip= new int[nmax];
  pointer= ip;
  for( int i= 0; i < nmax; ++i ) ip[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  char buf[33];
  sprintf( buf, "%i", nmax );
  leafs+= buf;
  leafs+= "]/I";
  brp= tp->Branch( label.c_str(), &ip[0], leafs.c_str(), 8000 );
  
}
void IntArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) { 

  const vector<int>* vp= (const vector<int>*) p;
  if( (int)vp->size() < nmax ) {
    std::cerr << "IntArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    for( int i= 0; i < nmax; ++i ) ((int*)pointer)[i]= (*vp)[i];
  }

}
IntDynArrCmsTreeColumn::IntDynArrCmsTreeColumn( const char* l, 
				  const vector<int> & v, 
				  const int & d, CmsTreeColumn* inp,
				  TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ), indexp( inp ) {

  // Make a new branch:
  int* np= (int*) indexp->getPointer();
  int* ip= new int[*np];
  pointer= ip;
  for( int i= 0; i < *np; ++i ) ip[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  leafs+= indexp->getLabel();
  leafs+= "]/I";
  brp= tp->Branch( label.c_str(), &ip[0], leafs.c_str(), 8000 );

}
void IntDynArrCmsTreeColumn::setDefValue() {

  if( pointer ) delete [] (int*)pointer;
  int nmax= *((int*)(indexp->getPointer()));
  int* ip= new int[nmax];
  pointer= ip;
  for( int i= 0; i < nmax; ++i ) ip[i]= defValue; 
  brp->SetAddress( &ip[0] );

}
void IntDynArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) {

  const vector<int>* vp= (const vector<int>*) p;
  int* np= (int*) cp->getPointer();
  if( *np > (int)vp->size() ) {
    std::cerr << "IntDynArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    int* ip= new int[*np];
    if( pointer ) delete[] (int*)pointer;
    pointer= ip;
    for( int i= 0; i < *np; ++i ) ip[i]= (*vp)[i];
    brp->SetAddress(  &ip[0] );
  }

}

// Float columns:
FloatCmsTreeColumn::FloatCmsTreeColumn( const char* l, const float & v, 
			  const float & d, 
			  TTree* tp ) : 
  CmsTreeColumn( l ), defValue( d ) {
				 
  // Create a new branch:
  pointer= new float;
  *(float*)pointer= v;
  std::string leafs( l ) ;
  leafs+= "/F";
  brp= tp->Branch( label.c_str(), pointer, leafs.c_str(), 8000 );

}
FloatArrCmsTreeColumn::FloatArrCmsTreeColumn( const char* l, 
				const vector<float> & v, 
				const float & d, TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ) {

  // Create a new branch:
  nmax= v.size();
  float* fp= new float[nmax];
  pointer= fp;
  for( int i= 0; i < nmax; ++i ) fp[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  char buf[33];
  sprintf( buf, "%i", nmax );
  leafs+= buf;
  leafs+= "]/F";
  brp= tp->Branch( label.c_str(), &fp[0], leafs.c_str(), 8000 );
  
}
void FloatArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) { 

  const vector<float>* vp= (const vector<float>*) p;
  if( (int)vp->size() < nmax ) {
    std::cerr << "FloatArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    for( int i= 0; i < nmax; ++i ) ((float*)pointer)[i]= (*vp)[i];
  }

}
FloatDynArrCmsTreeColumn::FloatDynArrCmsTreeColumn( const char* l, 
				      const vector<float> & v, 
				      const float & d, CmsTreeColumn* ip,
				      TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ), indexp( ip ) {

  // Make a new branch:
  int* np= (int*) indexp->getPointer();
  float* fp= new float[*np];
  pointer= fp;
  for( int i= 0; i < *np; ++i ) fp[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  leafs+= indexp->getLabel();
  leafs+= "]/F";
  brp= tp->Branch( label.c_str(), &fp[0], leafs.c_str(), 8000 );

}
void FloatDynArrCmsTreeColumn::setDefValue() {

  if( pointer ) delete [] (float*)pointer;
  int nmax= *((int*)(indexp->getPointer()));
  float* fp= new float[nmax];
  pointer= fp;
  for( int i= 0; i < nmax; ++i ) fp[i]= defValue; 
  brp->SetAddress( &fp[0] );

}
void FloatDynArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) {

  const vector<float>* vp= (const vector<float>*) p;
  int* np= (int*) cp->getPointer();
  if( *np > (int)vp->size() ) {
    std::cerr << "IntDynArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    float* fp= new float[*np];
    if( pointer ) delete[] (float*)pointer;
    pointer= fp;
    for( int i= 0; i < *np; ++i ) fp[i]= (*vp)[i];
    brp->SetAddress(  &fp[0] );
  }

}


// Double columns:
DoubleCmsTreeColumn::DoubleCmsTreeColumn( const char* l, const double & v, 
			    const double & d, TTree* tp ) : 
  CmsTreeColumn( l ), defValue( d ) {
				 
  // Create a new branch:
  pointer= new double;
  *(double*)pointer= v;
  std::string leafs( l ) ;
  leafs+= "/D";
  brp= tp->Branch( label.c_str(), pointer, leafs.c_str(), 8000 );

}
DoubleArrCmsTreeColumn::DoubleArrCmsTreeColumn( const char* l, 
				  const vector<double> & v, 
				  const double & d, TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ) {

  // Create a new branch:
  nmax= v.size();
  double* dp= new double[nmax];
  pointer= dp;
  for( int i= 0; i < nmax; ++i ) dp[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  char buf[33];
  sprintf( buf, "%i", nmax );
  leafs+= buf;
  leafs+= "]/D";
  brp= tp->Branch( label.c_str(), &dp[0], leafs.c_str(), 8000 );
  
}
void DoubleArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) { 

  const vector<double>* vp= (const vector<double>*) p;
  if( (int)vp->size() < nmax ) {
    std::cerr << "DoubleArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    for( int i= 0; i < nmax; ++i ) ((double*)pointer)[i]= (*vp)[i];
  }

}
DoubleDynArrCmsTreeColumn::DoubleDynArrCmsTreeColumn( const char* l, 
					const vector<double> & v, 
					const double & d, CmsTreeColumn* ip,
					TTree* tp ) :
  CmsTreeColumn( l ), defValue( d ), indexp( ip ) {

  // Make a new branch:
  int* np= (int*)(indexp->getPointer());
  double* dp= new double[*np];
  pointer= dp;
  for( int i= 0; i < *np; ++i ) dp[i]= v[i];
  std::string leafs( label.c_str() );
  leafs+= "[";
  leafs+= indexp->getLabel();
  leafs+= "]/D";
  brp= tp->Branch( label.c_str(), &dp[0], leafs.c_str(), 8000 );

}
void DoubleDynArrCmsTreeColumn::setDefValue() {

  if( pointer ) delete [] (double*)pointer;
  int nmax= *((int*)(indexp->getPointer()));
  double* dp= new double[nmax];
  pointer= dp;
  for( int i= 0; i < nmax; ++i ) dp[i]= defValue; 
  brp->SetAddress( &dp[0] );

}
void DoubleDynArrCmsTreeColumn::setValue( const void* p, CmsTreeColumn* cp ) {

  const vector<double>* vp= (const vector<double>*) p;
  int* np= (int*) cp->getPointer();
  if( *np > (int)vp->size() ) {
    std::cerr << "IntDynArrCmsTreeColumn::setValue: input vector too short, "
	 << "use default values" << std::endl;
    setDefValue();
  }
  else {
    double* dp= new double[*np];
    if( pointer ) delete[] (double*)pointer;
    pointer= dp;
    for( int i= 0; i < *np; ++i ) dp[i]= (*vp)[i];
    brp->SetAddress(  &dp[0] );
  }

}


