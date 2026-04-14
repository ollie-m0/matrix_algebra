// We aim to invert an n dimensional (square) matrix via Gaussian elimination.

#include <iostream>

void gaussian_elimination(double* arr_ptr, double* inv_arr_return, int n)
{

    // Define the identity matrix, which will be returned as the inverted matrix at the end.
    double inv_arr[n*n];
    inv_arr_return = inv_arr;
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            inv_arr[n*i+j] = (i == j) ? 1 : 0;
        }
    }

    // Reconstruct the initial array from the pointer.
    double arr[n*n];
    for (int i=0; i<n*n; i++)
    {
        arr[i] = *(arr_ptr + i);
    }



    // Reduce matrix to row-echelon form.
    double factor;
    for (int k=0; k<n-1; k++)
    {
        for (int i=k+1; i<n; i++)
        {
            factor = -(arr[i*n + k] / arr[k*(n+1)]); // Factor by which to multiply kth row before addition to ith row.
            if (factor==0) { continue; }  // While the float comparrison is potentially dangerous, this line is just to skip
                                          // code that doesn't do anything anyway, so not a big deal if the condition fails.
            // std::cout << "factor " << factor << std::endl;
            // std::cout << "arr[k..] " << arr[k*(n+1)] << std::endl;
            // std::cout << "arr[i..] " << arr[i*n + k] << std::endl;
            for (int j=0; j<n; j++)
            {
                arr[i*n + j]     = arr[i*n + j]     + factor * arr[k*n + j];
                inv_arr[i*n + j] = inv_arr[i*n + j] + factor * inv_arr[k*n + j];
            }
        }
    }

    // Back substitution
    for (int k=n-1; k>0; k--)
    {
        for (int i=0; i<k; i++)
        {
            factor = -(arr[i*n + k] / arr[k*(n+1)]);
            for (int j=0; j<n; j++)
            {
                arr[i*n + j]     = arr[i*n + j]     + factor * arr[k*n + j];
                inv_arr[i*n + j] = inv_arr[i*n + j] + factor * inv_arr[k*n + j];
            }
        }
    }

    // Divide by relevant factors such that the "old" matrix becomes the identity
    for (int i=0; i<n; i++)
    {
        factor = 1 / arr[n*i + i];
        for (int j=0; j<n; j++)
        {
            arr[n*i + j] = arr[n*i + j] * factor;
            inv_arr[n*i + j] = inv_arr[n*i + j] * factor;
        }
    }




    std::cout << "inverted matrix:" << std::endl;

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            std::cout << " " << inv_arr[i*n + j] << " ";
            if (j==n-1) { std::cout << std::endl; }
        }
    }


}
