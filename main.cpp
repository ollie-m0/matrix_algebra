#include <iostream>
#include "inversion.hpp"

int main()
{
    int n = 4;
    double* arr_ptr;
    double arr[n*n] = {  1,  3,  4,  3,
                        -1, -3,  7,  8,
                         5,  6,  6,  4,
                         1,  2,  2,  6
                       };
    arr_ptr = arr;

    double inverse[n*n];
    double* inverse_ptr;
    inverse_ptr = inverse;

    gaussian_elimination(arr, inverse_ptr, n);

    return 0;
}
