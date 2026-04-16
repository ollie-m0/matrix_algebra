#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>

void print_array(double* arr, int n, std::string title="", int setw_size=10)
{
    for (int i=0; i<(setw_size+1)*n; i++)
        std::cout << "_";
    std::cout << "__" << std::endl;
    std::cout << title << std::endl;

    for (int i=0; i<n; i++)
    {
        std::cout << '[';
        for (int j=0; j<n; j++)
        {
            std::cout << ' ' << std::setw(setw_size) << *(arr + i*n + j);
            if (j==n-1) { std::cout << ']' << std::endl; }
        }
    }
    for (int i=0; i<(setw_size+1)*n; i++)
        std::cout << "_";
    std::cout << "__" << std::endl;
}

void change_zero_row(double* arr, double* inv_arr, int n, int row, int column, double epsilon=1e-9)
{
    for (int j=row+1; j<n; j++)
    {
        if (std::fabs(*(arr + j*n + column)) <= epsilon)
        {
            if (j!=n-1) { continue; }
            else
            {
                std::cout << "Error encountered: divide by zero unavoidable. Check that det(matrix) != 0" << std::endl;
                std::abort();
            }
        }
        else
        {
            for (int k=0; k<n; k++)
            {
                *(arr + row*n + k)     += *(arr + j*n + k);
                *(inv_arr + row*n + k) += *(inv_arr + j*n + k);
            }
            break;
        }
    }
}



void gaussian_elimination(double* arr_ptr, double* inv_arr_return, int n, double epsilon=1e-9)
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

    // Handle the case in which the 0th element is 0
    if (std::fabs(arr[0]) <= epsilon)
    {
        change_zero_row(&arr[0], &inv_arr[0], n, 0, 0, epsilon);
    }


    // Reduce matrix to row-echelon form.
    double factor;
    for (int k=0; k<n-1; k++)
    {

        if (std::fabs(arr[k*(n+1)]) <= epsilon)
        // Avoid divide-by-zero errors
        {
            change_zero_row(&arr[0], &inv_arr[0], n, k, k, epsilon);
        }

        // Matrix reduction
        for (int i=k+1; i<n; i++)
        {
            factor = -(arr[i*n + k] / arr[k*(n+1)]); // Factor by which to multiply kth row before addition to ith row.

            if (factor==0.0) { continue; }  // While the float comparrison is potentially dangerous, this line is just to skip
                                            // code that wouldn't do anything anyway, so not a big deal if the condition fails.
            for (int j=0; j<n; j++)
            {
                arr[i*n + j]     += factor * arr[k*n + j];
                inv_arr[i*n + j] += factor * inv_arr[k*n + j];
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
            inv_arr[n*i + j] = inv_arr[n*i + j] * factor;
            // No need to perform the same opperation on arr[n*i + j],
            // as we are finished with that array.
        }
    }

    // Print the matrix (only for small n)
    if (n<=50)
    { print_array(&inv_arr[0], n, "Inverted matrix:", 10); }

}
