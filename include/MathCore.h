#ifndef _MATH_CORE_H
#define _MATH_CORE_H

#include "presicion.h"
//#include "blitz/array.h"

//Array<real,2> AxisXtoCell(const Array<real,2>&, const Array<real,2>&);
//Array<real,2> AxisYtoCell(const Array<real,2>&, const Array<real,2>& );


void insertionSort(int a[], int b[], int c[], real d[], int array_size)
{
     int i, j, index, index2, index3;
     real index4;
     
     for (i = 1; i < array_size; ++i)
     {
          index  = a[i];
          index2 = b[i];
          index3 = c[i];
          index4 = d[i];
          for (j = i; j > 0 && a[j-1] > index; j--){
               a[j] = a[j-1];
               b[j] = b[j-1];
               c[j] = c[j-1];
               d[j] = d[j-1];
          }
          a[j] = index;
          b[j] = index2;
          c[j] = index3;
          d[j] = index4;
     }
}

void add_ptr(int a[], int sizeA, int b[], int sizeB)
{

    for(int i = 0 ; i < sizeB ; i++){
        a[b[i]]++;
    } 


    for(int i = 1 ; i < sizeA ; i++){
        a[i] += a[i-1];
        //sum = a[i];
    } 

    for(int i = sizeA-1; i >= 0 ; i--){
         a[i+1] = a[i];
        //sum = a[i];
    } 
    
    a[0] = 0; // start value 



}

#endif
