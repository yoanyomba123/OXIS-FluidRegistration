#ifndef __ImageContiguousArray_h
#define __ImageContiguousArray_h

#include <iostream>
#include <iomanip>
#include <Constants.hpp>
#include <stdlib.h>

using namespace std;

/*
* ImageContiguousArray.h
*
* Class provides an in memory data structure with various initialization 
* schemes, operations, and bounds checking functionalities
*
* Data in this data struture is stored in a contiguous block thus, operations 
* on its data members must be handled in a specialized manner
*
* Given an input image of size NxN, this data structures allocated a singular
* contiguous block of memory with the size being N (rows) x N (cols). Element 
* access is specified by indicating an initial offset which in this case
* signifies the primal element in a row and an index delimiting the offset 
* from the first element in that row.
*
* <pre>
*   - data access : object(row_size x i + col_index)
* </pre>
*
* @Author: D Yoan L Mekontchou Yomba
* @Date: 01/13/2019
*/

class ImageContiguousArray{
  // define data members
  protected: 
    double numRows; // number of rows
    double numCols; // number of cols
    double totalSize; // complete size of image
    double* rows; // rows pointer
    double* cols; // cols pointer
    long *data; // array data
  
  // define constructor
  ImageContiguousArray(){
    numRows = Constants::GridConstants::rows;
    numCols = Constants::GridConstants::cols;
    totalSize = numRows * numCols;
    rows = &numRows;
    cols = &numCols;
    // allocate a default amount of memory
    allocateMemoryAndZeroOut(numRows, numCols);

  }

  ImageContiguousArray(double rows, double cols){
    numRows = rows;
    numCols = cols;
    totalSize = numRows * numCols;
    rows = &numRows;
    cols = &numCols;
    allocateMemoryAndZeroOut(rows, cols);
  }


  ImageContiguousArray(Image2D& image){
      numRows = image->rows;
      numCols = image->cols;
      totalSize = numRows * numCols;
      rows = &numRows;
      cols = &numCols;
      allocateMemoryAndInitiliaze(rows, cols, image);
  }

  virtual ~ImageContiguousArray(){
    freeMemory();
  }


  private:
    void allocateMemoryAndZeroOut(double rows, double cols){
        data = (long *) calloc(rows * cols , sizeof(long));
        if(data == NULL){
          fprintf(stderr, "Fatal: failed to allocate %zu bytes. \n", (rows * cols * sizeof(long))); 
        }
    }

    void allocateMemoryAndInitiliaze(double rows, double cols,Image2D& image){
      data = (long *) calloc(rows * cols , sizeof(long));
      if(data == NULL){
        fprintf(stderr, "Fatal: failed to allocate %zu bytes. \n", (rows * cols * sizeof(long))); 
      }

      int count = 1;
      int index = 0;
      long imageData = image->imagePixelGrid;

      // allocate values from 2D array to 1D memory block
      for(int i =0; i < totalSize; i++){
        if(count % rows == 0){
          index += 1;
          count = 1
        }
        *(data+i) = imageData[index][count - 1]
        count += 1;
      }
    }

    void freeMemory(){
      free(data);
    }
  
  public:

    // returns the maximal pixel value of the image
    long getMaximalPixelValue(){
      long maximalValue = 0;
      for(int i = 0; i < totalSize; i++){
        if(maximalValue > *(data+i)) 
          maximalValue = *(data+i);
      }
      return maximalValue;
    }

    // returns the minimal pixel value of the image
    long getMinimalPixelValue(){
      long minimalValue;
      for(int i = 0; i < totalSize; i++){
        
        if(i == 0 )
          minimalValue = *(data+i);
        
        if(minimalValue > *(data+i)) 
          minimalValue = *(data+i);
      }
      return minimalValue;
    }


    // gets the pixel at a given location
    long getPixelValue(int row, int col){
      return *(array * row + col);
    }


}

#endif