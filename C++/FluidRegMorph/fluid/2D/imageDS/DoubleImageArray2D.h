#ifndef __DoubleImageArray2D__
#define __DoubleImageArray2D__

#include <iostream>
#include <iomanip>
using namespace std;

/* DoubleImageArray2D.h
 *
 * Class provides a two dimensional data structure with initialization schemes
 * numeric operations and bounds checking.
 *
 * <pre>
 * - The primal index is by default 0 in accordance to the typical C convention
 * - Data for the array is stored in rows
 * - Element access is performed by using the following construct Object(i,j) where (i,j) signify the ith row and the jth column in the object data structure
 * </pre>
 *
 * For each instance of the class, a copy constructor is utilized to create a duplicate instance in order to protect inproper or accidental deletion
 *
 * @Author D Yoan L Mekontchou Yomba (Drexel U.)
 * @Date 01/13/2019
 */

class DoubleImageArray2D {

    protected:
        double* dataPtr; // data pointer
        long idx1Begin;  // starting index for the initial coordinate (rows)
        long idx2Begin;  // starting index for the secondary coordinate (cols)
        long idx1End;    // ending index for the initial coordiante (rows)
        long idx2End;    // ending index for the secondary coordiante (cols)
        long idx1Size;   // primary coordinate size (rows)
        long idx2Size;   // secondary coordiante size (cols)


    // constructor definition
    DoubleImageArray2D(){
        dataPtr = 0;
        idx1Size = 0;
        idx2Size = 0;
        idx1Begin = 0;
        idx2Begin = 0;
        idx1End = 0;
        idx2End = 0;
    }


    // constructor
    DoubleImageArray2D(long m, long n){
        dataPtr = 0;
        idx1Size = 0;
        idx2Size = 0;
        idx1Begin = 0;
        idx2Begin = 0;
        idx1End = 0;
        idx2End = 0;
        initialize(m,n);
    }

    // constructor
    DoubleImageArray2D(double *d, long m, long n){
        dataPtr = 0;
        idx1Size = 0;
        idx2Size = 0;
        idx1Begin = 0;
        idx2Begin = 0;
        idx1End = 0;
        idx2End = 0;
        initialize(d,m,n);
    }

    // constructor
    DoubleImageArray2D(const DoubleImageArray2D &D){

        idx1Size = D.idx1Size;
        idx2Size = D.idx2Size;
        idx1Begin = D.idx1Begin;
        idx2Begin = D.idx2Begin;
        idx1End = D.idx1End;
        idx2End = D.idx2End;

        dataPtr = new double[idx1Size * idx2Size];
        for(long i; i < idx1Size * idx2Size; i++){
            dataPtr[i] = D.dataPtr[i];
        }
    }

    // class destructor
    virtual ~DoubleImageArray2D(){
        if(dataPtr != 0){
            delete [] dataPtr;
        }
    }

    /* 
    * def - initiliaze the data structure more specifically initialize the data pointer to 0
    *   input:
    *       @param m - the number of rows
    *       @param n - the number of cols
    */
    void initialize(long m, long n){
        if(idx1Size != m || idx2Size != n){
            delete [] dataPtr;
            dataPtr = new double[m * n];
        }
        idx1Size = m;
        idx2Size = n;
        idx1Begin = 0;
        idx2Begin = 0;
        idx1End = idx1Begin + (idx1Size - 1);
        idx2End = idx2Begin + (idx2Size - 1);

        for (long i =0; i < idx1Size * idx2Size; i++){
            dataPtr[i] = 0.0;
        }
    }


    /* 
    * def -  initialize the data structure as well as the data pointer based on the passed in parameter ref
    *   input:
    *       @param d - the data pointer to initialize this structure to
    *       @param m - the number of rows
    *       @param n -  the number of cols
    */
    void initialize(double* d, long m, long n){
        initialize(m, n);
        for(long i =0; i < idx1Size * idx2Size; i++){
            dataPtr[i] = d[i];
        }
    }

    /* 
    * def - gets the dataPointer of this structure
    *       output:
    *           @return double* - double Pointer
    */
    double* getDataPointer(){
        return dataPtr;
    }

    /*
    *   def - set the start index 
    *      input:
    *           @param i - the index to be set
    */
    void setIdx1Begin(long i){
        idx1Begin = i;
        idx1End = idx1Begin + (idx1Size - 1);
    }

    /*
    *   def - set the start index 
    *      input:
    *           @param i - the index to be set
    */
    void setIdx2Begin(long i){
        idx2Begin = i;
        idx2End = idx2Begin + (idx2Size - 1);
    }

    // return the begin index for the data structurce (idx1 == row)
    long getIdx1Begin() const {
        return idx1Begin;
    }

    // returns the begin index for the data structure (idx2 == col)
    long getIdx2Begin() const {
        return idx2Begin;
    }

    // returns the magnitude of rows in the data structure (idx1 == row) 
    long getIdx1Size() const {
        return idx1Size;
    }

    // returns the magnitude of cols in the data structure (idx2 == col)
    long getIdx2Size() const {
        return idx2Size;
    }

    // returns the ending index of the rows in the data structure (idx1 == rows)
    long getIdx1End() const {
        return idx1End;
    }

    // returns the ending index of the cols in the data structure (idx2 == cols)
    long getIdx2End() const {
        return idx2End;
    }

    /*
    *   def - sets all elements in the data pointer to the input param
    *       input:
    *           @param d -  the value to set the data pointer to
    */
    void setDataPointerToValue(double d){
        for (long i = 0; i < idx1Size * idx2Size; i ++){
            dataPtr[i] = d;
        }
    }

    /*
    *   def - takes the dot product of the input and the current class instance
    *       input: 
    *           @param &D - the memory addres of the input array structure
    *       ouput:
    *           @result result - the long value obtained by taking the dot product of the two matrixes
    */
    double dotProduct(const DoubleImageArray2D &D){
        long result = 0.0;
        for(long i = 0; i < idx1Size * idx2Size; i++){
            result += D.dataPtr[i] * dataPtr[i];
        }
        return result;
    }

    // override the default c++ operator and multiplies this instance's dataPtr by a multiplier
    void operator*=(double multiplier){
        for(long i = 0; i < idx1Size * idx2Size; i++){
            dataPtr[i] *= multiplier;
        }
    }

    // negate the input arra's dataPtr from this instance's
    void operator-=(const DoubleImageArray2D &D){

        // check if this class instance's coordinates (1 or 2) are zero initialize to that of the input array
        if(idx1Size*idx2Size == 0){
            initialize(D.idx1Size, D.idx2Size);
        }
        
        for(long i =0; i < idx1Size*idx2Size; i++){
            dataPtr[i] -= D.dataPtr[i];
        }
    }

    // add the data pointer value from the input array to this classes' instance
    void operator+=(const DoubleImageArray2D &D){
        // check if this class instance's coordinates (1 or 2) are zero initialize to that of the input array
        if(idx1Size*idx2Size == 0){
            initialize(D.idx1Size, D.idx2Size);
        }
        
        for(long i =0; i < idx1Size*idx2Size; i++){
            dataPtr[i] += D.dataPtr[i];
        }
    }

    // check for equality between two arays
    void operator=(const DoubleImageArray2D &D){
        if(idx1Size*idx2Size == 0){
            initialize(D.idx1Size, D.idx2Size);
            idx1Size = D.idx1Size;
            idx2Size = D.idx2Size;
            idx1Begin = D.idx1Begin;
            idx2Begin = D.idx2Begin;
            idx1End = D.idx1End;
            idx2End = D.idx2End;
        }
        for(long i = 0; i < idx1Size*idx2Size; i++){
            dataPtr[i] = D.dataPtr[i];
        }
    }

    // perfoms division of all elements in a dataPtr by an input value
    DoubleImageArray2D operator/(double value){
        DoubleImageArray2D R(*this);
        for(long i =0; i < idx1Size*idx2Size; i++){
            R.dataPtr[i] /= value;
        }
        return R;
    }

    // multiplies this instance with the input value provided as an input param
    DoubleImageArray2D operator*(double value){
        DoubleImageArray2D R(*this);
        for(long i =0; i , idx1Size*idx2Size; i++){
            R.dataPtr[i] *= value;
        }
        return R;
    }

    // mulitplies data pointer values of the input array by the 1st param which is a double
    friend DoubleImageArray2D operator*(double value, const DoubleImageArray2D &D){
        DoubleImageArray2D R(D);
        for(long i  = 0; i < idx1Size * idx2Size; i++){
            R.dataPtr[i] *= value;
        }
        return R;
    }

    // overload the negation operator to negate a value from each dataPtr element
    DoubleImageArray2D operator-(double value){
        DoubleImageArray2D R(*this);
        for(long i =0; i < idx1Size*idx2Size; i++){
            R.dataPtr[i] -= value;
        }
        return R;
    }

    // overload the negation operator to negate from one array by another input array
    DoubleImageArray2D operator-(const DoubleImageArray2D &D){
        DoubleImageArray2D R(*this);
        for(long i =0; i < idx1Size*idx2Size; i++){
            R.dataPtr[i] -= D.dataPtr[i];
        }
        return R;
    }

    // overload the plus operator to add a given value to all elements in the dataPtr
    DoubleImageArray2D operator+(double value){
        DoubleImageArray2D R(*this);
        for(long i =0; i < idx1Size*idx2Size; i++){
            R.dataPtr[i] += value;
        }
        return R;
    }

    // overload the plus operator for matrix addition
    DoubleImageArray2D operator+(const DoubleImageArray2D &D){
        DoubleImageArray2D R(*this);
        for(long i =0; i < idx1Size*idx2Size; i++){
            R.dataPtr[i] += D.dataPtr[i];
        }
        return R;
    }

    friend ostream&  operator <<(ostream& outStream, const DoubleImageArray2D& A){
        long i; long j;
        for(j = A.index2End; j >=  A.index2Begin; j--){
            for(i = A.index1Begin; i <=  A.index1End; i++){
                outStream <<  setw(5) << A(i,j) << " ";
            }
            outStream << endl;
        }
        return outStream; 
    }
};

#endif