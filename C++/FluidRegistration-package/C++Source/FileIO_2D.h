#ifndef __FileIO_2D__
#define __FileIO_2D__

#include <vector>
#include "DoubleArray1D.h"
#include "DoubleArray2D.h"
#include "Grid2D.h"

typedef vector<    int    > vector1Dint;
typedef vector<vector1Dint> vector2Dint;
typedef vector<vector2Dint> vector3Dint;
typedef vector<vector3Dint> vector4Dint;
typedef vector<vector4Dint> vector5Dint;

typedef vector<  double   > vector1Ddouble;

//
//#####################################################################
//						FileIO_2D.h
//#####################################################################
//
// Igor Yanovsky (C) UCLA
// version:  11/25/2008
//
//#####################################################################
//


class FileIO_2D
{
public:

	// PGM ASCII Info:
	void getImageInfo( string filename, long& m, long& n );

	// PGM ASCII READ:
	DoubleArray2D readPGM( string filename );
//	DoubleArray2D FileIO_2D::readPGM( string filename, long& m, long& n, const long bx, const long by );
	DoubleArray2D readPGM( string filename, long& m, long& n, const long bx, const long by );

	// PGM ASCII WRITE:
	void writePGM(   const DoubleArray2D& A, const string& str );
	void writePGM(   const DoubleArray2D& A, const string& str, const long& stepCount );
	void writePGM(   const DoubleArray2D& A, const string& str, const long& stepCount, const long bx, const long by );
	void writePGM(   const DoubleArray2D& A, const string& str, const long& stepCount, const long& b );
	void writePGMsc( const DoubleArray2D& A, const string& str );
	void writePGMsc( const DoubleArray2D& AA, const string& str, const long& stepCount );
	void writePGMsc( const DoubleArray2D& AA, const string& str, const long& stepCount, const long bx, const long by);
	void writePGMsc( const DoubleArray2D& AA, long r, long l, long t, long b, const string& str );
	void writePGMsc( const DoubleArray2D& AA, const string& str, long r, long l, long t, long b, const long& stepCount );

	// PPM ASCII READ:
	// not implemented yet:
	void readPPM ( DoubleArray2D& A, DoubleArray2D& B, DoubleArray2D& C, string filename );

	// ASCII READ:
	void readDAT2D( DoubleArray2D& A, long m, long n, string filename );
	void readDAT2D( DoubleArray2D& A, string filename );
	void readDAT1D( DoubleArray1D& A, string filename );

	// ASCII WRITE:
	void write_ascii( const vector2Dint&   A, const string& str );
	void write_ascii( const DoubleArray2D& A, const string& str, const long& stepCount );
	void write_ascii( const DoubleArray2D& A, const string& str, const long& stepCount, const long& w );
	void write_ascii( const DoubleArray2D& A, const string& str );
	void write_ascii( const DoubleArray2D& A, const string& str, const long& stepCount, const long bx, const long by );
	void write_ascii( const DoubleArray1D& A, const string& str, const long& stepCount );
	void write_ascii( const DoubleArray1D& A, const string& str );
	void write_ascii( const vector1Ddouble& A, const string& str );

	// BINARY READ:
	void read_bin_uc(  string filename, DoubleArray2D& A );		// unsigned char
	void read_bin_usi( string filename, DoubleArray2D& A, bool byteswap );
	void read_bin_ssi( string filename, DoubleArray2D& A, bool byteswap );
	void read_bin_float( string filename, DoubleArray2D& A );

	// BINARY WRITE:
	void write_bin_usi( const DoubleArray2D& A, const string& str, long factor ); // unsigned short int
	void write_bin_usi( const DoubleArray2D& A, const string& str, const long& stepCount, long factor );
	void write_bin_usi( const DoubleArray2D& A, const string& str, const long& stepCount, long factor, const long bx, const long by );

	void write_bin_ssi( const DoubleArray2D& A, const string& str, long factor );
	void write_bin_float( const DoubleArray2D& A, const string& str );
	
	// HEADER:
	void make_header( const DoubleArray2D& A, const string& str );
	void make_header( const DoubleArray2D& A, const string& str, const long& bx, const long& by );
};

#endif

