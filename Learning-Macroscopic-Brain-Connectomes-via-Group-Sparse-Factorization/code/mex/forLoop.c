#include "mex.h"
#include <string.h>
#include <math.h>
#include <time.h>

// THESE are bad names, should be changed. They indicate the max number of active a and f for a voxel
#define MAXA 1000
#define MAXF 100

// These are the total number of possible atoms, fasiclces and voxels
#define TOTALA 200000
#define TOTALFAS 1100
#define TOTALV 20000

// Indexing array used for all voxels; must be initialized to zero to start
int NA[TOTALA];
int NF[TOTALFAS];
int NV[TOTALV];

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

// Column major ordering in Matlab (unlike c), so index is col * num_rows + row
// For a matrix of num_rows x num_cols, return the index for (row, col)
unsigned long get_index(unsigned long row, unsigned long col, unsigned long num_rows, unsigned long num_cols) {
      return (col * num_rows + row);
}     

// sign function
double sign(double a) {
    double threshold = 0.000001;
    if(a > threshold) return 1;
    if(a < -threshold) return -1;
    return 0;
}

// return the max possible voxel index in Gv
// the index is the voxel index in original Phi
unsigned long get_max_Gv(double * sub, unsigned long row, unsigned long col) {
    unsigned long max = 0;
    unsigned long i;
    double temp;
    for(i = 0; i < row; ++i) {
       temp = sub[get_index(i, 0, row, col)];
       if(temp > max) {
           max = temp;
       }
    }
    return (unsigned long)max;
}

// store an sptensor into cell-array according to dimension dim
// this dimension dim corresponds to the Nvb dimension
// e.g. sub(some_row, dim) = a, where a \in [1, Nvb]
// The only exception is for Gv, where we are on the first dim of Nv
// but when we process in the future, we only take out the Nv corresponding to 
// Nvb, with the help of voxels
// the returned cell array has element like {sub, val, [row], [col]}
// sub stores the positions of nonzero elements, 
// val stores the correspoinding nonzero values
// [row] and [col] are dimension of the sparse tensor
// O(# of nonzero elements)
void fill_cell_array(mxArray * array, unsigned long Nvb, double * sub, double * val, unsigned long row, unsigned long col, unsigned long dim) {
    unsigned long i;
    unsigned long * count = mxCalloc(Nvb, sizeof(unsigned long));
    mxArray * sub_array;
    mxArray * rc_pair_one;
    double * rc_pt_one;
    mxArray * rc_pair_two;
    double * rc_pt_two;
    // 4 elements in each array for each voxel, which is the element of the final array
    mwSize sub_length = 4;
    
    // calculate the count of each voxels of Nvb appears in sub
    // namely nonzero in the sparse tensor
    for(i = 0; i < row; ++i) {
      ++count[(unsigned long)(sub[get_index(i, dim, row, col)] - 1)]; 
    }
    
    double * sub_temp;
    double * val_temp;
    unsigned long my_row;
    unsigned long my_col;
    // since it is a slice of voxels of the original tensor, so the column reduces 1
    my_col = col - 1;
    
    for(i = 0; i < Nvb; ++i) {
        sub_array = mxCreateCellArray(1, &sub_length);
        rc_pair_one = mxCreateDoubleMatrix(1, 1, mxREAL); 
        rc_pt_one = mxGetPr(rc_pair_one);
        rc_pair_two = mxCreateDoubleMatrix(1, 1, mxREAL); 
        rc_pt_two = mxGetPr(rc_pair_two);
        my_row = count[i];
        // if no such voxel
        if(my_row < 1) {
            //mxSetCell(sub_array, 0, 0);
            //mxSetCell(sub_array, 1, 0);
            //rc_pt_one[0] = 0;
            //rc_pt_two[0] = 0;
            mxSetCell(sub_array, 0, mxCreateDoubleMatrix(1, my_col, mxREAL));
            mxSetCell(sub_array, 1, mxCreateDoubleMatrix(1, 1, mxREAL));
            rc_pt_one[0] = 1;
            rc_pt_two[0] = my_col;
            
            mxSetCell(sub_array, 2, rc_pair_one);
            mxSetCell(sub_array, 3, rc_pair_two);
            
            sub_temp = mxGetPr(mxGetCell(sub_array, 0));
            val_temp = mxGetPr(mxGetCell(sub_array, 1));
            for (unsigned long j = 0; j < my_col; ++j) {
                sub_temp [get_index(0, j, 1, my_col)] = 1;
            }
            val_temp[get_index(0, 0, 1, 1)] = 0;
            
            mxSetCell(array, i, sub_array);
            continue;
        }
        // set the dimensions of the returned objects
        mxSetCell(sub_array, 0, mxCreateDoubleMatrix(my_row, my_col, mxREAL));
        mxSetCell(sub_array, 1, mxCreateDoubleMatrix(my_row, 1, mxREAL));
        rc_pt_one[0] = my_row;
        rc_pt_two[0] = my_col;
        mxSetCell(sub_array, 2, rc_pair_one);
        mxSetCell(sub_array, 3, rc_pair_two);
        mxSetCell(array, i, sub_array);      
    }
    
    unsigned long n;
    unsigned long j, k;
    unsigned long * index = mxCalloc(Nvb, sizeof(unsigned long));
    // iterate the row of sub (each nonzero element's index)
    // set the values of the returned objects
    for(i = 0; i < row; ++i) {
       // n is the index of the voxel from dimension dim
       n = (unsigned long)(sub[get_index(i, dim, row, col)] - 1);
       sub_array = mxGetCell(array, n);
       sub_temp = mxGetPr(mxGetCell(sub_array, 0));
       val_temp = mxGetPr(mxGetCell(sub_array, 1));
       my_row = count[n];
       k = 0;
       for(j = 0; j < col; ++j) {
          if(j == dim) {
              continue; // since dim is the position of voxel, the result ignores it
          }
          sub_temp[get_index(index[n], k, my_row, my_col)] = sub[get_index(i, j, row, col)]; 
          ++k;
       }
       val_temp[index[n]] = val[i];
       ++index[n];
    }
    mxFree(count);
    mxFree(index);
}

int check_active_one(double index, double *indices, int num_ind) {
    // For loop over each array and determine valid indices
    int i;
    for (i = 0; i < num_ind; i++) {
        if (indices[i] == index) {
            return 1;
        }
    }
    // If reached end, index_one not a valid index, so index cannot be active
    return 0;
}


int check_active(double index_one, double index_two, double * atomind, double *fasind, int num_atomind, int num_fasind) {
     // Faster to loop atomind, so do that first (because orientations are smaller than fascicles ?)
    if (check_active_one(index_one, atomind, num_atomind) == 0) {
        return 0;
    }
    return check_active_one(index_two, fasind, num_fasind);
}

// hadmard product for sparse tensor
// e.g. a .* b
// O(#(a) * #(b))
void h_product_spt(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, char flag) {
    unsigned long i;
    unsigned long j;
    unsigned long m;
    unsigned long k;
    unsigned long s;
    unsigned long n;
    *res_row = 0;
    *res_col = col_one;
    for(i = 0; i < row_one; ++i) {       
        for(j = 0; j < row_two; ++j) {           
            for(k = 0; k < col_two; ++k) {
                if(sub_one[get_index(i, k, row_one, col_one)] != sub_two[get_index(j, k, row_two, col_two)]) {
                    break;
                }
            }
            if(k == col_two) {
                (*res_row)++;
            }                    
        }
    }
    m = 0;
    for(i = 0; i < row_one; ++i) {       
        for(j = 0; j < row_two; ++j) {           
            for(k = 0; k < col_two; ++k) {
                if(sub_one[get_index(i, k, row_one, col_one)] != sub_two[get_index(j, k, row_two, col_two)]) {
                    break;
                }
            }
            if(k == col_two) {
                for(s = 0; s < col_one; ++s) {
                    res_sub[get_index(m, s, *res_row, *res_col)] = sub_one[get_index(i, s, row_one, col_one)];  
                }
                if (flag == 0) {
                    res_val[m] = val_one[i] * val_two[j];
                }
                else if(flag == 1) {
                    res_val[m] = val_one[i] * sign(val_two[j]);
                }
                ++m;  
            }                    
        }
    }
}

// this is for outer product of two sparse tensor
// col_one = 1, col_two = 1 (both tensors have only one dimension)
// O(#(a) * #(b))
void ttt_outer_product_spt(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two) {
    unsigned long i;
    unsigned long j;
    unsigned long m;
    unsigned long row_count = 0;
   
    *res_row = row_one * row_two;
    *res_col = col_one + col_two;
   
    m = 0;
   
    for(i = 0; i < row_one; ++i) {       
        for(j = 0; j < row_two; ++j) {
            //dim_two can be only 0 or 1
            res_sub[get_index(m, 0, *res_row, *res_col)] = sub_one[get_index(i, 0, row_one, col_one)];
            res_sub[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(j, 0, row_two, col_two)];
            res_val[m] = val_one[i] * val_two[j];
            ++m;
        }
    }
}

// sparse tensor multiplies another sparse tensor (absolute value)
// e.g. ttt(A,fabs(B(:,v,:)),1,1)
// O(#(a) * #(b))
void ttt_one(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two) {
    unsigned long i;
    unsigned long j;
    unsigned long m;
    unsigned long k;
    unsigned long s;
    unsigned long n;
    double temp;
    *res_row = 0;
    *res_col = col_one + col_two - 2;
    for(i = 0; i < row_one; ++i) {       
        temp = sub_one[get_index(i, 0, row_one, col_one)];
        for(j = 0; j < row_two; ++j) {            
            if(sub_two[get_index(j, 0, row_two, col_two)] == temp) {
                (*res_row)++;
            }                  
        }
    }
   
    double * res_sub_t = mxCalloc(((*res_row) * (*res_col)), sizeof(double));
    double * res_val_t = mxCalloc((*res_row), sizeof(double));
    
    m = 0;
    for(i = 0; i < row_one; ++i) {      
        temp = sub_one[get_index(i, 0, row_one, col_one)];
        for(j = 0; j < row_two; ++j) {            
            if(sub_two[get_index(j, 0, row_two, col_two)] == temp) {
                res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_one[get_index(i, 1, row_one, col_one)];
                res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(j, 1, row_two, col_two)];
                res_val_t[m] = val_one[i] * fabs(val_two[j]);
                ++m;
            }                  
        }
    }
    double first, second;
    m = 0;
    for(i = 0; i < (*res_row); ++i) {
        first = res_sub_t[get_index(i, 0, *res_row, *res_col)];
        second = res_sub_t[get_index(i, 1, *res_row, *res_col)];
        for(j = 0; j < m; ++j) {
            if(first == res_sub_t[get_index(j, 0, *res_row, *res_col)] && second == res_sub_t[get_index(j, 1, *res_row, *res_col)]) {
                res_val_t[j] += res_val_t[i];
                break;
            }
        }
        if(j == m) {
            res_sub_t[get_index(m, 0, *res_row, *res_col)] = first;
            res_sub_t[get_index(m, 1, *res_row, *res_col)] = second;
            res_val_t[m] = res_val_t[i];
            ++m;
        }
    }
   
   //res_sub = (double *)malloc(sizeof(double)*(m * (*res_col)));   
   //res_val = (double *)malloc(sizeof(double)*m);
   
    for(i = 0; i < m; ++i) {
        res_sub[get_index(i, 0, m, *res_col)] = res_sub_t[get_index(i, 0, *res_row, *res_col)];
        res_sub[get_index(i, 1, m, *res_col)] = res_sub_t[get_index(i, 1, *res_row, *res_col)];
        res_val[i] = res_val_t[i];
    }
    *res_row = m;
    mxFree(res_sub_t);
    mxFree(res_val_t);
}

/*
// sparse tensor multiplies another sparse tensor
// e.g. ttt(A(v,:),B,1,2)
void ttt_two(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, double v)
{
   unsigned long i;
   unsigned long j;
   unsigned long m;
   unsigned long k;
   unsigned long s;
   unsigned long n;
   double temp;
   *res_row = 0;
   *res_col = col_one+col_two-3;
    for(i = 0; i<row_one; ++i)
    {       
        if(sub_one[get_index(i,0,row_one,col_one)]==v)
        {
            temp = sub_one[get_index(i,1,row_one,col_one)];
            for(j=0;j<row_two;++j)
            {
                if(sub_two[get_index(j,1,row_two,col_two)]==temp)
                {
                    (*res_row)++;
                } 
            }
        }
    }
       
    double * res_sub_t = calloc(((*res_row) * (*res_col)),sizeof(double)); 
    double * res_val_t = calloc((*res_row), sizeof(double));
   
    m = 0;
    for(i = 0; i<row_one; ++i)
    {       
        if(sub_one[get_index(i,0,row_one,col_one)]==v)
        {
            temp = sub_one[get_index(i,1,row_one,col_one)];
            for(j=0;j<row_two;++j)
            {
                if(sub_two[get_index(j,1,row_two,col_two)]==temp)
                {
                    res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_two[get_index(j,0,row_two,col_two)];
                    res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(j,2,row_two,col_two)];
                    res_val_t[m] = val_one[i] * val_two[j];
                    ++m;
                } 
            }
        }
    }
    
    double first, second;
    m = 0;
    for(i=0;i<(*res_row);++i)
    {
       first = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       second = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       for(j=0;j<m;++j)
       {
        if(first == res_sub_t[get_index(j, 0, *res_row, *res_col)] && second == res_sub_t[get_index(j, 1, *res_row, *res_col)])
        {
            res_val_t[j] += res_val_t[i];
            break;
        }
       }
       if(j==m)
       {
         res_sub_t[get_index(m, 0, *res_row, *res_col)] = first;
         res_sub_t[get_index(m, 1, *res_row, *res_col)] = second;
         res_val_t[m] = res_val_t[i];
         ++m;
       }
    }
   
    //res_sub = (double *)malloc(sizeof(double)*(m * (*res_col)));   
    //res_val = (double *)malloc(sizeof(double)*m);
   
   for(i = 0;i<m;++i)
   {
       res_sub[get_index(i, 0, m, *res_col)] = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       res_sub[get_index(i, 1, m, *res_col)] = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       res_val[i] = res_val_t[i];
   }
   
   *res_row = m;
   free(res_sub_t);
   free(res_val_t);
}
*/

// sparse tensor multiplies another sparse tensor
// the first tensor is one dimensional, the second one is three dimensional
// e.g. ttt(A,B,1,2)
// additionally takes a subset of atomind and fasind, to restrict the atom and fasin in B
void ttt_two(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, double *fasind, int num_fasind, int flag_usescreening) {
    unsigned long i;
    unsigned long j;
    unsigned long m;
    unsigned long k;
    unsigned long s;
    unsigned long n;
    double temp;
    *res_row = 0;
    *res_col = col_one + col_two - 2;
    for(j = 0; j < row_two; ++j) {
       // Ignore any indices in B that are not in fasind and atomind
        if(flag_usescreening && check_active_one(sub_two[get_index(j,2,row_two,col_two)], fasind, num_fasind) == 0) {
            continue;
        }
      
        for(i = 0; i < row_one; ++i) {
            if(sub_two[get_index(j, 1, row_two, col_two)] == sub_one[get_index(i, 0, row_one, col_one)]) {
                (*res_row)++;
            } 
        }
    }
   // MARTHA: Old way, changed to above to allow subsets of fasind and atomind
   /*for(i = 0; i<row_one; ++i)
    {       
        temp = sub_one[get_index(i,0,row_one,col_one)];
        for(j=0;j<row_two;++j)
        {
            if(sub_two[get_index(j,1,row_two,col_two)]==temp)
            {
                (*res_row)++;
            } 
        }
        }*/
       
    double * res_sub_t = calloc(((*res_row) * (*res_col)),sizeof(double)); 
    double * res_val_t = calloc((*res_row), sizeof(double));

    /* Which of these two ways is faster? I think its method 1*/
    /*
    int method = 1;
    m = 0;
    if (method == 0) {
    for(i = 0; i<row_one; ++i)
    {       
        temp = sub_one[get_index(i,0,row_one,col_one)];
        for(j=0;j<row_two;++j)
        {
           // Only other change for subset is to add a check_active here
            if(sub_two[get_index(j,1,row_two,col_two)]==temp
               && (!flag_usescreening || check_active_one(sub_two[get_index(j,2,row_two,col_two)], fasind, num_fasind)))
            {
                res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_two[get_index(j,0,row_two,col_two)];
                res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(j,2,row_two,col_two)];
                res_val_t[m] = val_one[i] * val_two[j];
                ++m;
            } 
        }
    }
    }*/
    //if (method == 1) {
    for(j = 0; j < row_two; ++j) {
      // Ignore any indices in B that are not in fasind and atomind
        if (flag_usescreening && check_active_one(sub_two[get_index(j, 2, row_two, col_two)], fasind, num_fasind) == 0) {
            continue;
        }
        for(i = 0; i < row_one; ++i) {
           // Only other change for subset is to add a check_active here
            if (sub_two[get_index(j, 1, row_two, col_two)] == sub_one[get_index(i, 0, row_one, col_one)]) {
                res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_two[get_index(j, 0, row_two, col_two)];
                res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(j, 2, row_two, col_two)];
                res_val_t[m] = val_one[i] * val_two[j];
                ++m;
            }
        }
    }
    
    double first, second;
    m = 0;
    for(i = 0; i < (*res_row); ++i) {
       first = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       second = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       for(j = 0; j < m; ++j) {
           if(first == res_sub_t[get_index(j, 0, *res_row, *res_col)] && second == res_sub_t[get_index(j, 1, *res_row, *res_col)]) {
               res_val_t[j] += res_val_t[i];
               break;
           }
       }
       if(j == m) {
         res_sub_t[get_index(m, 0, *res_row, *res_col)] = first;
         res_sub_t[get_index(m, 1, *res_row, *res_col)] = second;
         res_val_t[m] = res_val_t[i];
         ++m;
       }
    }
   
   for(i = 0;i < m; ++i) {
       res_sub[get_index(i, 0, m, *res_col)] = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       res_sub[get_index(i, 1, m, *res_col)] = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       res_val[i] = res_val_t[i];
   }
   
   *res_row = m;
   free(res_sub_t);
   free(res_val_t);
}

// sparse tensor multiplies another sparse tensor
// the first tensor is one dimensional, the second one is three dimensional
// e.g. ttt(A,B,1,2)
// additionally takes a subset of atomind and fasind, to restrict the atom and fasin in B
// extra bookkeeping to avoid for looping
// Assume a is an indicator vector
void ttt_two_fast(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, double *fasind, int num_fasind, int flag_usescreening) {
   unsigned long i;
   unsigned long j;
   unsigned long m;
   *res_row = 0;
   *res_col = col_one + col_two - 2;

   double A[MAXA][MAXF];
   unsigned long INDICESA[MAXA];
   unsigned long INDICESFAS[MAXF];
   memset(A, 0, sizeof(double) * MAXA * MAXF);
  
       
   unsigned long fas_ind_b, voxel_ind_b, voxel_ind_a, or_ind_b;
   int current_or_index_A = 1, or_loc;
   int current_fas_index_A = 1, fas_loc;

   // Assume NA, NF and Nv were zeroed last time to save computation
   // Fill in Nv 
    for(i = 0; i < row_one; ++i){
        voxel_ind_a = sub_one[get_index(i, 0, row_one, col_one)];
        NV[voxel_ind_a] = 1;
    }

    // Fill in indexing for fascicles
    // Not making this a global, since TOTALFAS is small
    int ACTIVE_F[TOTALFAS];
    memset(ACTIVE_F, 0, sizeof(int) * TOTALFAS);
    for(i = 0; i < num_fasind; ++i){
        ACTIVE_F[(int)fasind[i]] = 1;
    }
      
    for(j = 0; j < row_two; ++j) {
        voxel_ind_b = sub_two[get_index(j, 1, row_two, col_two)];
        fas_ind_b = sub_two[get_index(j, 2, row_two, col_two)];

        // if (NV[voxel_ind_b] < 1 || check_active_one(fas_ind_b, fasind, num_fasind) == 0)
        if (NV[voxel_ind_b] < 1 || ACTIVE_F[fas_ind_b] < 1) { 
            continue;
        }

        // Otherwise, this is a valid index in B
      
        // Add to array; if already active, simply sum up in A
        // Otherwise, add a new index
        or_ind_b = sub_two[get_index(j, 0, row_two, col_two)];
        or_loc = current_or_index_A;
        fas_loc = current_fas_index_A;
        if (NA[or_ind_b] > 0) { 
            or_loc = NA[or_ind_b];
        }
        else {
            NA[or_ind_b] = current_or_index_A;
            INDICESA[current_or_index_A] = or_ind_b;
            current_or_index_A++;
        }

        if (NF[fas_ind_b] > 0) {
            fas_loc = NF[fas_ind_b];
        }
        else {
            NF[fas_ind_b] = current_fas_index_A;
            INDICESFAS[current_fas_index_A] = fas_ind_b;                     
            current_fas_index_A++;
        }
        // Assume val_one[i] = 1
        // A[or_loc][fas_loc] += val_one[i] * val_two[j];
        A[or_loc][fas_loc] += val_two[j];
    }
    m = 0;
    for(i = 1; i < current_or_index_A; ++i) {
       for(j = 1; j < current_fas_index_A; ++j) {
          // If A[i][j] is zero then that index was not active
          if (A[i][j] > 0) {
             ++m;
          }
       }
    }

    *res_row = m;
    m = 0;
    for(i = 1;i < current_or_index_A; ++i) {
       for(j = 1;j < current_fas_index_A; ++j) {
          // If A[i][j] is zero then that idex was not active
          if (A[i][j] > 0) {
             res_sub[get_index(m, 0, *res_row, *res_col)] = INDICESA[i];
             res_sub[get_index(m, 1, *res_row, *res_col)] = INDICESFAS[j];
             res_val[m] = A[i][j];
             ++m;
          }
       }
    }
    // Zero out NA and NF  and Nv
    for(i = 1; i < current_or_index_A; ++i) {
        NA[INDICESA[i]] = 0;
    }
    for(j = 1; j < current_fas_index_A; ++j) {
        NF[INDICESFAS[j]] = 0;
    }
    for(i = 0; i < row_one; ++i) {
     voxel_ind_a = sub_one[get_index(i, 0, row_one, col_one)];
     NV[voxel_ind_a] = 0;
    }    
}

// sparse tensor multiplies another sparse tensor (absolute value), 
// only for two dimensions, namely dim_one and dim_two \in {0,1} 
// e.g. ttt(A,B,2,1)
void ttt(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, unsigned long dim_one, unsigned long dim_two) {
   unsigned long i;
   unsigned long j;
   unsigned long m;
   unsigned long k;
   unsigned long s;
   unsigned long n;
   double temp;
   *res_row = 0;
   *res_col = col_one + col_two - 2;
   for(i = 0; i < row_one; ++i) {
       temp = sub_one[get_index(i, dim_one, row_one, col_one)];
       for(j = 0; j < row_two; ++j) {
           if(sub_two[get_index(j, dim_two, row_two, col_two)] == temp) {
               (*res_row)++;
           }
       }
   }
   //double * res_sub_t = (double *)malloc(sizeof(double)*((*res_row) * (*res_col)));
   //double * res_val_t = (double *)malloc(sizeof(double)*(*res_row));
   double * res_sub_t = calloc(((*res_row) * (*res_col)), sizeof(double)); 
   double * res_val_t = calloc((*res_row), sizeof(double));
   m = 0;
   for(i = 0; i < row_one; ++i) {
       temp = sub_one[get_index(i, dim_one, row_one, col_one)];
       for(j = 0; j < row_two; ++j) {
           if(sub_two[get_index(j, dim_two, row_two, col_two)] == temp) {
               res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_one[get_index(i, 1 - dim_one, row_one, col_one)];
               res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(j, 1 - dim_two, row_two, col_two)];
               res_val_t[m] = val_one[i] * fabs(val_two[j]);
               ++m;
           }
       }
   }
   double first, second;
   m = 0;
   for(i = 0; i < (*res_row); ++i) {
       first = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       second = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       for(j = 0; j < m; ++j) {
           if(first == res_sub_t[get_index(j, 0, *res_row, *res_col)] && second == res_sub_t[get_index(j, 1, *res_row, *res_col)]) {
               res_val_t[j] += res_val_t[i];
               break;
           }
       }
       if(j == m) {
         res_sub_t[get_index(m, 0, *res_row, *res_col)] = first;
         res_sub_t[get_index(m, 1, *res_row, *res_col)] = second;
         res_val_t[m] = res_val_t[i];
         ++m;
       }
    }
   //res_sub = (double *)malloc(sizeof(double)*(m * (*res_col)));   
   //res_val = (double *)malloc(sizeof(double)*m);
   for(i = 0; i < m; ++i) {
       res_sub[get_index(i, 0, m, *res_col)] = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       res_sub[get_index(i, 1, m, *res_col)] = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       res_val[i] = res_val_t[i];
   }
   *res_row = m;
   free(res_sub_t);
   free(res_val_t);
}


// add three sparse tensor
// e.g. A(:,v,:) = A(:,v,:)+B+lambda*C
// but additionally restrict the atoms and fascicles considered, where ATOM_IND=0 and FAS_IND = 2

void add_subset(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, double * sub_three, double * val_three, unsigned long row_three, unsigned long col_three, double lambda, double * atomind, double *fasind, int num_atomind, int num_fasind) {

    unsigned long i;
    unsigned long j;
    unsigned long m = 0;
    
    *res_row = row_one + row_two + row_three;
    *res_col = col_one;
    double * res_sub_t = calloc((*res_row) * (*res_col), sizeof(double));
    double * res_val_t = calloc(*res_row, sizeof(double));
    
    for(i = 0; i < row_one; ++i) {
       // Check if index active
        if (check_active(sub_one[get_index(i, 0, row_one, col_one)], sub_one[get_index(i, 1, row_one, col_one)], atomind, fasind, num_atomind, num_fasind)) {
            res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_one[get_index(i, 0, row_one, col_one)];
            res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_one[get_index(i, 1, row_one, col_one)];
            res_val_t[m] = val_one[i];
            ++m;
        }
    }
    
    for(i = 0; i < row_two; ++i) {
        // If index not active, not need to consider it further
        if (check_active(sub_two[get_index(i, 0, row_two, col_two)], sub_two[get_index(i, 1, row_two, col_two)], atomind, fasind, num_atomind, num_fasind) == 0) {
             continue;
        }

        // Else, its a valid index and we can add it to the tensor
        for(j = 0; j < m; ++j) {
            if(res_sub_t[get_index(j, 0, *res_row, *res_col)] == sub_two[get_index(i, 0, row_two, col_two)] && res_sub_t[get_index(j, 1, *res_row, *res_col)] == sub_two[get_index(i, 1, row_two, col_two)]) {
                res_val_t[j] += val_two[i];
                break;
            }
        }
        
        if(j == m) {
            res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_two[get_index(i, 0, row_two, col_two)];
            res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(i, 1, row_two, col_two)];
            res_val_t[m] = val_two[i];
            ++m;
        }
    }
    
    for(i = 0; i < row_three; ++i) {
       // If index not active, no need to consider it further
        if (check_active(sub_three[get_index(i, 0, row_three, col_three)], sub_three[get_index(i, 1, row_three, col_three)], atomind, fasind, num_atomind, num_fasind) == 0) {
            continue;
        }
        for(j = 0; j < m; ++j) {
            if(res_sub_t[get_index(j, 0, *res_row, *res_col)] == sub_three[get_index(i, 0, row_three, col_three)] && res_sub_t[get_index(j, 1, *res_row, *res_col)] == sub_three[get_index(i, 1, row_three, col_three)]) {
                res_val_t[j] += lambda * val_three[i];
                //mexPrintf("Entereted!: %lf, %lf, %lf\n", res_val_t[j], lambda, val_three[i]);
                break;
            }
        }
        if(j == m) {
            res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_three[get_index(i, 0, row_three, col_three)];
            res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_three[get_index(i, 1, row_three, col_three)];
            res_val_t[m] = lambda * val_three[i];
            ++m;
        }
    }
    for(i = 0; i < m; ++i) {
        res_sub[get_index(i, 0, m, *res_col)] = res_sub_t[get_index(i, 0, *res_row, *res_col)];
        res_sub[get_index(i, 1, m, *res_col)] = res_sub_t[get_index(i, 1, *res_row, *res_col)];
        res_val[i] = res_val_t[i];
        //mexPrintf("Outer loop!: %lf\n", res_val[i]);
    }
    *res_row = m;
    free(res_sub_t);
    free(res_val_t);
}

// add three sparse tensors 
// e.g. A(:,v,:) = A(:,v,:) + B + lambda * C
void add(double * res_sub, double * res_val, unsigned long * res_row, unsigned long * res_col, double * sub_one, double * val_one, unsigned long row_one, unsigned long col_one, double * sub_two, double * val_two, unsigned long row_two, unsigned long col_two, double * sub_three, double * val_three, unsigned long row_three, unsigned long col_three, double lambda) {
    unsigned long i;
    unsigned long j;
    unsigned long m;
    
    *res_row = row_one + row_two + row_three;
    *res_col = col_one;
    //double * res_sub_t = (double *)malloc(sizeof(double)*((*res_row) * (*res_col)));
    //double * res_val_t = (double *)malloc(sizeof(double)*(*res_row));
    double * res_sub_t = calloc((*res_row) * (*res_col), sizeof(double));
    double * res_val_t = calloc(*res_row, sizeof(double));
    
    m = row_one;
    for(i = 0; i < m; ++i) {
        res_sub_t[get_index(i, 0, *res_row, *res_col)] = sub_one[get_index(i, 0, row_one, col_one)];
        res_sub_t[get_index(i, 1, *res_row, *res_col)] = sub_one[get_index(i, 1, row_one, col_one)];
        res_val_t[i] = val_one[i];
    }
    
    for(i = 0; i < row_two; ++i) {
        for(j = 0; j < m; ++j) {
            if(res_sub_t[get_index(j ,0, *res_row, *res_col)] == sub_two[get_index(i, 0, row_two, col_two)] && res_sub_t[get_index(j, 1, *res_row, *res_col)] == sub_two[get_index(i, 1, row_two, col_two)]) { 
                res_val_t[j] += val_two[i];
                break;
            }
        }
        
        if(j == m) {
            res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_two[get_index(i, 0, row_two, col_two)];
            res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_two[get_index(i, 1, row_two, col_two)];
            res_val_t[m] = val_two[i];
            ++m;
        }
    }
    
    for(i = 0; i < row_three; ++i) {
        for(j = 0; j < m; ++j) {
           if(res_sub_t[get_index(j, 0, *res_row, *res_col)] == sub_three[get_index(i, 0, row_three, col_three)] && res_sub_t[get_index(j, 1, *res_row, *res_col)] == sub_three[get_index(i, 1, row_three, col_three)]) {
               res_val_t[j] += lambda * val_three[i];
               break;
           }
        }
        if(j == m) {
            res_sub_t[get_index(m, 0, *res_row, *res_col)] = sub_three[get_index(i ,0, row_three, col_three)];
            res_sub_t[get_index(m, 1, *res_row, *res_col)] = sub_three[get_index(i, 1, row_three, col_three)];
            res_val_t[m] = lambda * val_three[i];
            ++m;
        }
    }
    
   for(i = 0; i < m; ++i) {
       res_sub[get_index(i, 0, m, *res_col)] = res_sub_t[get_index(i, 0, *res_row, *res_col)];
       res_sub[get_index(i, 1, m, *res_col)] = res_sub_t[get_index(i, 1, *res_row, *res_col)];
       res_val[i] = res_val_t[i];
   }
    
   *res_row = m;
   free(res_sub_t);
   free(res_val_t);
}

void print_array(double * subs, double * vals, int numrows, int numcols) {
    int i;
    for (i = 0; i < numrows; i++) {
        mexPrintf("%lf %lf: %lf\n", subs[get_index(i, 0, numrows, numcols)], subs[get_index(i, 1, numrows, numcols)], vals[i]);
    }
    mexPrintf("End\n");
}

// The inputs are lambda_g1,a_mask, f_mask, entry_masks, Ga, Gv, 
// grad_p1_x_t (namely B'(Ydiff)), grad_g1_x2 (namely 1/A), Phib_ls, sign_Phib_ls (do not pass it, we can calculate it from Phib_ls)
// voxels (index for voxel, this is a matrix (Nv*1) rather than a sparse tensor)
// g (already contains some part of gradient and we need add result to that)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // If optimize fascicles, then allows fascicles to become active if 
    // they are active in the neighbors of the voxels and lets objective reduce this;
    // this case correpsonds to f_mask;
    // Else, fascicles are fixed, and new ones cannot become active;
    // this case corresponds to f_screen
    int flag_optimizefascicles = 0;

    // Initialize NA and NF and Nv
    memset(NA, 0, sizeof(NA[0]) * TOTALA);
    memset(NF, 0, sizeof(NF[0]) * TOTALFAS);
    memset(NV, 0, sizeof(NV[0]) * TOTALV);

    double lambda_g1 = mxGetScalar(prhs[0]);
   
    double * sub_aMask = mxGetPr(prhs[1]);
    double * val_aMask = mxGetPr(prhs[2]);
    double * sub_fMask = mxGetPr(prhs[3]);
    double * val_fMask = mxGetPr(prhs[4]);
    double * sub_entryMask = mxGetPr(prhs[5]);
    double * val_entryMask = mxGetPr(prhs[6]);
    double * sub_Ga = mxGetPr(prhs[7]);
    double * val_Ga = mxGetPr(prhs[8]);  
    double * sub_Gv = mxGetPr(prhs[9]);
    double * val_Gv = mxGetPr(prhs[10]);

    double * sub_grad_p1_xt = mxGetPr(prhs[11]);
    double * val_grad_p1_xt = mxGetPr(prhs[12]);
    double * sub_A_inv = mxGetPr(prhs[13]);
    double * val_A_inv = mxGetPr(prhs[14]);
    double * sub_Phib = mxGetPr(prhs[15]);
    double * val_Phib = mxGetPr(prhs[16]);
    double * voxels = mxGetPr(prhs[17]);
    double * sub_g = mxGetPr(prhs[18]);
    double * val_g = mxGetPr(prhs[19]);  

    int SCREEN_IND = 20;
    double * sub_aScreen = mxGetPr(prhs[SCREEN_IND]);
    double * val_aScreen = mxGetPr(prhs[SCREEN_IND + 1]);
    double * sub_fScreen = mxGetPr(prhs[SCREEN_IND + 2]);
    double * val_fScreen = mxGetPr(prhs[SCREEN_IND + 3]);
   
    unsigned long Nvb = mxGetM(prhs[17]); // voxels contains Nvb indexes for voxels
    mwSize length = Nvb;
   
    unsigned long grad_p1_xt_row = mxGetM(prhs[11]);
    unsigned long grad_p1_xt_col = mxGetN(prhs[11]);
    mxArray * grad_p1_xt_array = mxCreateCellArray(1, &length);
   
    unsigned long aMask_row = mxGetM(prhs[1]);
    unsigned long aMask_col = mxGetN(prhs[1]);
    mxArray * aMask_array = mxCreateCellArray(1, &length);
   
    unsigned long fMask_row = mxGetM(prhs[3]);
    unsigned long fMask_col = mxGetN(prhs[3]);
    mxArray * fMask_array = mxCreateCellArray(1, &length);
   
    unsigned long entryMask_row = mxGetM(prhs[5]);
    unsigned long entryMask_col = mxGetN(prhs[5]);
    mxArray * entryMask_array = mxCreateCellArray(1, &length);
   
    unsigned long Ga_row = mxGetM(prhs[7]);
    unsigned long Ga_col = mxGetN(prhs[7]);
   
    unsigned long Phib_row = mxGetM(prhs[15]);
    unsigned long Phib_col = mxGetN(prhs[15]);
    mxArray * Phib_array = mxCreateCellArray(1, &length);
   
    unsigned long Gv_row = mxGetM(prhs[9]);
    unsigned long Gv_col = mxGetN(prhs[9]);
   
    mwSize Gv_length;
    Gv_length = get_max_Gv(sub_Gv, Gv_row, Gv_col);;
    mxArray * Gv_array = mxCreateCellArray(1, &Gv_length);
   
    unsigned long A_inv_row = mxGetM(prhs[13]);
    unsigned long A_inv_col = mxGetN(prhs[13]);
    
   
    unsigned long g_row = mxGetM(prhs[18]);
    unsigned long g_col = mxGetN(prhs[18]); 
    mxArray * g_array = mxCreateCellArray(1, &length);

    unsigned long aScreen_row = mxGetM(prhs[SCREEN_IND]);
    unsigned long aScreen_col = mxGetN(prhs[SCREEN_IND]);
    mxArray * aScreen_array = mxCreateCellArray(1, &length);
   
    unsigned long fScreen_row = mxGetM(prhs[SCREEN_IND + 2]);
    unsigned long fScreen_col = mxGetN(prhs[SCREEN_IND + 2]);
    mxArray * fScreen_array = mxCreateCellArray(1, &length);
   
    unsigned long MAX_ROWS = 100000;
    unsigned long MAX_COLS = 100;
   
    double *sub_grad_p1_xv;
    double *val_grad_p1_xv;  
    unsigned long grad_p1_xv_row = MAX_ROWS;
    unsigned long grad_p1_xv_col = MAX_COLS;
    sub_grad_p1_xv = mxCalloc(grad_p1_xv_row * grad_p1_xv_col, sizeof(double));
    val_grad_p1_xv = mxCalloc(grad_p1_xv_row, sizeof(double));

    double * sub_grad_p1_v;
    double * val_grad_p1_v;
    unsigned long grad_p1_v_row = MAX_ROWS;
    unsigned long grad_p1_v_col = MAX_COLS;
    sub_grad_p1_v = mxCalloc(grad_p1_v_row * grad_p1_v_col, sizeof(double));
    val_grad_p1_v = mxCalloc(grad_p1_v_row, sizeof(double));
   
    double * sub_grad_p1_v_new;
    double * val_grad_p1_v_new;
    unsigned long grad_p1_v_new_row = MAX_ROWS;
    unsigned long grad_p1_v_new_col = MAX_COLS;
    sub_grad_p1_v_new = mxCalloc(grad_p1_v_new_row * grad_p1_v_new_col, sizeof(double));
    val_grad_p1_v_new = mxCalloc(grad_p1_v_new_row, sizeof(double));
   
    double * sub_grad_g1_x;
    double * val_grad_g1_x;
    unsigned long grad_g1_x_row = MAX_ROWS;
    unsigned long grad_g1_x_col = MAX_COLS;
    sub_grad_g1_x = mxCalloc(grad_g1_x_row * grad_g1_x_col, sizeof(double));
    val_grad_g1_x = mxCalloc(grad_g1_x_row, sizeof(double));
   
    double * sub_grad_g1_x3;
    double * val_grad_g1_x3;
    unsigned long grad_g1_x3_row = MAX_ROWS;
    unsigned long grad_g1_x3_col = MAX_COLS;
    sub_grad_g1_x3 = mxCalloc(grad_g1_x3_row * grad_g1_x3_col, sizeof(double));
    val_grad_g1_x3 = mxCalloc(grad_g1_x3_row, sizeof(double));
   
    double * sub_grad_g1_x4;
    double * val_grad_g1_x4;
    unsigned long grad_g1_x4_row = MAX_ROWS;
    unsigned long grad_g1_x4_col = MAX_COLS;
    sub_grad_g1_x4 = mxCalloc(grad_g1_x4_row * grad_g1_x4_col, sizeof(double));
    val_grad_g1_x4 = mxCalloc(grad_g1_x4_row, sizeof(double));
   
    double * sub_temp;
    double * val_temp;
    unsigned long temp_row = MAX_ROWS;
    unsigned long temp_col = MAX_COLS;
    sub_temp = mxCalloc(temp_row * temp_col, sizeof(double));
    val_temp = mxCalloc(temp_row, sizeof(double));
   
    double * sub_grad_g1_v;
    double * val_grad_g1_v;
    unsigned long grad_g1_v_row = MAX_ROWS;
    unsigned long grad_g1_v_col = MAX_COLS;
    sub_grad_g1_v = mxCalloc(grad_g1_v_row * grad_g1_v_col, sizeof(double));
    val_grad_g1_v = mxCalloc(grad_g1_v_row, sizeof(double));
   
    double * sub_eachG;
    double * val_eachG;
    unsigned long default_eachG_row = MAX_ROWS;
    unsigned long default_eachG_col = MAX_COLS;
    unsigned long eachG_row = 0; // not used right away
    unsigned long eachG_col = 0;
    //sub_eachG = mxCalloc(eachG_row*eachG_col, sizeof(double));
    //val_eachG = mxCalloc(eachG_row, sizeof(double));
   
    double * sub_tempG;
    double * val_tempG;
    unsigned long tempG_row, tempG_col;
   
    double * sub_finalG;
    double * val_finalG;
    unsigned long finalG_row, finalG_col;
    
    unsigned long finalM = 0; // how many rows will the returned g have. 
   
    unsigned long vi, i;
    double v;

    unsigned long max_computed_size = 0;
   
    fill_cell_array(grad_p1_xt_array, Nvb, sub_grad_p1_xt, val_grad_p1_xt, grad_p1_xt_row, grad_p1_xt_col, 1);
    fill_cell_array(aMask_array, Nvb, sub_aMask, val_aMask, aMask_row, aMask_col, 1);
    fill_cell_array(entryMask_array, Nvb, sub_entryMask, val_entryMask, entryMask_row, entryMask_col, 1);   
    fill_cell_array(fMask_array, Nvb, sub_fMask, val_fMask, fMask_row, fMask_col, 1);
    fill_cell_array(Phib_array, Nvb, sub_Phib, val_Phib, Phib_row, Phib_col, 1);
    fill_cell_array(g_array, Nvb, sub_g, val_g, g_row, g_col, 1);
    fill_cell_array(Gv_array, Gv_length, sub_Gv, val_Gv, Gv_row, Gv_col, 0);
    fill_cell_array(aScreen_array, Nvb, sub_aScreen, val_aScreen, aScreen_row, aScreen_col, 1);
    fill_cell_array(fScreen_array, Nvb, sub_fScreen, val_fScreen, fScreen_row, fScreen_col, 1);
   
    double * sub_grad_p1_xt_v;
    double * val_grad_p1_xt_v;
    unsigned long grad_p1_xt_v_row;
    unsigned long grad_p1_xt_v_col;

    double * sub_aMask_v;
    double * val_aMask_v;
    unsigned long aMask_v_row;
    unsigned long aMask_v_col;
   
    double * sub_fMask_v;
    double * val_fMask_v;
    unsigned long fMask_v_row;
    unsigned long fMask_v_col;
   
    double * sub_entryMask_v;
    double * val_entryMask_v;
    unsigned long entryMask_v_row;
    unsigned long entryMask_v_col;

    double * sub_fScreen_v;
    double * val_fScreen_v;
    unsigned long fScreen_v_row;
    unsigned long fScreen_v_col;

    double * sub_aScreen_v;
    double * val_aScreen_v;
    unsigned long aScreen_v_row;
    unsigned long aScreen_v_col;
   
    double * sub_Phib_v;
    double * val_Phib_v;
    unsigned long Phib_v_row;
    unsigned long Phib_v_col;
   
    double * sub_g_v;
    double * val_g_v;
    unsigned long g_v_row;
    unsigned long g_v_col;
   
    double * sub_Gv_v;
    double * val_Gv_v;
    unsigned long Gv_v_row;
    unsigned long Gv_v_col;
   
    mxArray * EachG_array = mxCreateCellArray(1, &length);
   
    mxArray * rc_pair_one;
    double * rc_pt_one;
    mxArray * rc_pair_two;
    double * rc_pt_two;
   
    mxArray * sub_for_g_v;
    mxArray * val_for_g_v;
   
    mwSize sub_length = 4; 
   
    mxArray * sub_array;

    int max_computed_sizes[10];
    double time_spent[10];
    for (i = 0; i < 9; i++) {
        max_computed_sizes[i] = 0;
        time_spent[i] = 0;
    }

    clock_t begin;
    clock_t end;

   for(vi = 0; vi < Nvb; ++vi) {
       v = voxels[vi];// get the index of the voxels

       sub_array = mxGetCell(grad_p1_xt_array, vi);
       sub_grad_p1_xt_v = mxGetPr(mxGetCell(sub_array, 0));
       val_grad_p1_xt_v = mxGetPr(mxGetCell(sub_array, 1));
       grad_p1_xt_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       grad_p1_xt_v_col = mxGetScalar(mxGetCell(sub_array, 3));
       
       sub_array = mxGetCell(aMask_array, vi);
       sub_aMask_v = mxGetPr(mxGetCell(sub_array, 0));
       val_aMask_v = mxGetPr(mxGetCell(sub_array, 1));
       aMask_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       aMask_v_col = mxGetScalar(mxGetCell(sub_array, 3));
       
       sub_array = mxGetCell(fMask_array, vi);
       sub_fMask_v = mxGetPr(mxGetCell(sub_array, 0));
       val_fMask_v = mxGetPr(mxGetCell(sub_array, 1));
       fMask_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       fMask_v_col = mxGetScalar(mxGetCell(sub_array, 3));
       
       sub_array = mxGetCell(entryMask_array, vi);
       sub_entryMask_v = mxGetPr(mxGetCell(sub_array, 0));
       val_entryMask_v = mxGetPr(mxGetCell(sub_array, 1));
       entryMask_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       entryMask_v_col = mxGetScalar(mxGetCell(sub_array, 3));

       sub_array = mxGetCell(aScreen_array, vi);
       sub_aScreen_v = mxGetPr(mxGetCell(sub_array, 0));
       val_aScreen_v = mxGetPr(mxGetCell(sub_array, 1));
       aScreen_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       aScreen_v_col = mxGetScalar(mxGetCell(sub_array, 3));

       sub_array = mxGetCell(fScreen_array, vi);
       sub_fScreen_v = mxGetPr(mxGetCell(sub_array, 0));
       val_fScreen_v = mxGetPr(mxGetCell(sub_array, 1));
       fScreen_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       fScreen_v_col = 1;
      
       sub_array = mxGetCell(Phib_array, vi);
       sub_Phib_v = mxGetPr(mxGetCell(sub_array, 0));
       val_Phib_v = mxGetPr(mxGetCell(sub_array, 1));
       Phib_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       Phib_v_col = mxGetScalar(mxGetCell(sub_array, 3));
     
       sub_array = mxGetCell(g_array, vi);
       sub_g_v = mxGetPr(mxGetCell(sub_array, 0));
       val_g_v = mxGetPr(mxGetCell(sub_array, 1));
       g_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       g_v_col = mxGetScalar(mxGetCell(sub_array, 3));
      
       sub_array = mxGetCell(Gv_array, v - 1);
       sub_Gv_v = mxGetPr(mxGetCell(sub_array, 0));
       val_Gv_v = mxGetPr(mxGetCell(sub_array, 1));
       Gv_v_row = mxGetScalar(mxGetCell(sub_array, 2));
       Gv_v_col = mxGetScalar(mxGetCell(sub_array, 3));

       /************************ First Compute gradient of regularizer, which doesnt make new fascicles active ******************/

       // this is grad_g1_x = ttt(Ga, Phib_ls(:,vi,:) .* sign_Phib_ls(:,vi,:), 1, 1);
       // grad_g1_x1 = \|Phib_ls\| \times_1 \G_A^\top$
       begin = clock();
       ttt_one(sub_grad_g1_x, val_grad_g1_x, &grad_g1_x_row, &grad_g1_x_col, sub_Ga, val_Ga, Ga_row, Ga_col, sub_Phib_v, val_Phib_v, Phib_v_row, Phib_v_col);
       max_computed_size = max(grad_g1_x_row, max_computed_size);
       max_computed_sizes[3] = max(grad_g1_x_row, max_computed_sizes[3]);
       end = clock();
       time_spent[3] += (double)(end - begin) / CLOCKS_PER_SEC;
       //mexPrintf("time3, sizes: %d, %d\n", grad_g1_x_row, grad_g1_x_col);
       
       // this is grad_g1_x3 = ttt(Gv(v,:), grad_g1_x2, 1, 2);
       // Can only use screening here for fascicles
       // ttt_two(sub_grad_g1_x3, val_grad_g1_x3, &grad_g1_x3_row, &grad_g1_x3_col, sub_Gv, val_Gv, Gv_row, Gv_col, sub_A_inv, val_A_inv, A_inv_row, A_inv_col, v);
       begin = clock();
       // ttt_two(sub_grad_g1_x3, val_grad_g1_x3, &grad_g1_x3_row, &grad_g1_x3_col, sub_Gv_v, val_Gv_v, Gv_v_row, Gv_v_col, sub_A_inv, val_A_inv, A_inv_row, A_inv_col, sub_fScreen_v, fScreen_v_row, 1);
       ttt_two_fast(sub_grad_g1_x3, val_grad_g1_x3, &grad_g1_x3_row, &grad_g1_x3_col, sub_Gv_v, val_Gv_v, Gv_v_row, Gv_v_col, sub_A_inv, val_A_inv, A_inv_row, A_inv_col, sub_fScreen_v, fScreen_v_row, 1);        
       max_computed_size = max(grad_g1_x3_row, max_computed_size);
       max_computed_sizes[4] = max(grad_g1_x3_row, max_computed_sizes[4]);
       end = clock();
       time_spent[4] += (double)(end - begin) / CLOCKS_PER_SEC;
       //mexPrintf("time4, sizes2: %d, %d\n", grad_g1_x3_row, grad_g1_x3_col);
       /*
       if (vi == 0) {
        mexPrintf("time %d:\n", 4);
        for (int i = 0; i < grad_g1_x3_row; ++i) {
            mexPrintf("%lf, %lf, %lf\n", sub_grad_g1_x3[get_index(i, 0, grad_g1_x3_row, grad_g1_x3_col)], sub_grad_g1_x3[get_index(i, 1, grad_g1_x3_row, grad_g1_x3_col)], val_grad_g1_x3[get_index(i, 0, grad_g1_x3_row, 1)]);
        }
       }
        **/

       // this is grad_g1_x4 = grad_g1_x .* grad_g1_x3;
       begin = clock();
       h_product_spt(sub_grad_g1_x4, val_grad_g1_x4, &grad_g1_x4_row, &grad_g1_x4_col, sub_grad_g1_x, val_grad_g1_x, grad_g1_x_row, grad_g1_x_col, sub_grad_g1_x3, val_grad_g1_x3, grad_g1_x3_row, grad_g1_x3_col,0);
       max_computed_size = max(grad_g1_x4_row, max_computed_size);
       max_computed_sizes[5] = max(grad_g1_x4_row, max_computed_sizes[5]);
       end = clock();
       time_spent[5] += (double)(end - begin) / CLOCKS_PER_SEC;
       //mexPrintf("time5, sizes3: %d, %d\n", grad_g1_x4_row, grad_g1_x4_col);
       // these are grad_g1_v = ttt(Ga, grad_g1_x4, 2, 1) .* sign_Phib_ls(:,vi,:);
       // ttt(Ga, grad_g1_x4, 2, 1) part
       begin = clock();
       ttt(sub_temp, val_temp, &temp_row, &temp_col, sub_Ga, val_Ga, Ga_row, Ga_col, sub_grad_g1_x4, val_grad_g1_x4, grad_g1_x4_row, grad_g1_x4_col, 1, 0);
       max_computed_size = max(temp_row, max_computed_size);
       max_computed_sizes[6] = max(temp_row, max_computed_sizes[6]);
       end = clock();
       time_spent[6] += (double)(end - begin) / CLOCKS_PER_SEC;
       //mexPrintf("time6, sizes4: %d, %d\n", temp_row, temp_col);
       // .* sign_Phib_ls(:,vi,:) part
       begin = clock();
       h_product_spt(sub_grad_g1_v, val_grad_g1_v, &grad_g1_v_row, &grad_g1_v_col, sub_temp, val_temp, temp_row, temp_col, sub_Phib_v, val_Phib_v, Phib_v_row, Phib_v_col,1);  
       max_computed_size = max(grad_g1_v_row, max_computed_size);  
       max_computed_sizes[7] = max(grad_g1_v_row, max_computed_sizes[7]);  
       end = clock();
       time_spent[7] += (double)(end - begin) / CLOCKS_PER_SEC;
       //mexPrintf("time7, sizes: %d, %d\n", grad_g1_v_row, grad_g1_v_col);
       /************************ Second Compute gradient of loss, which can make new fascicles active ******************/
       
       // this is grad_p1_x_v = grad_p1_x_t(:,vi) .* a_mask(:,vi);
       // MARTHA: this can be further reduced to only the atoms for a voxel, rather than atom vincinity
       // so I have done so. We are not changing the oreitnations, so the loss update would ignore those
       begin = clock();
       
       //h_product_spt(sub_grad_p1_xv, val_grad_p1_xv, &grad_p1_xv_row, &grad_p1_xv_col, sub_grad_p1_xt_v, val_grad_p1_xt_v, grad_p1_xt_v_row, grad_p1_xt_v_col, sub_aMask_v, val_aMask_v, aMask_v_row, aMask_v_col, 0);
       h_product_spt(sub_grad_p1_xv, val_grad_p1_xv, &grad_p1_xv_row, &grad_p1_xv_col, sub_grad_p1_xt_v, val_grad_p1_xt_v, grad_p1_xt_v_row, grad_p1_xt_v_col, sub_aScreen_v, val_aScreen_v, aScreen_v_row, aScreen_v_col, 0);       
       
       max_computed_size = max(grad_p1_xv_row, max_computed_size);
       max_computed_sizes[0] = max(grad_p1_xv_row, max_computed_sizes[0]);
       end = clock();
       time_spent[0] += (double)(end - begin) / CLOCKS_PER_SEC;
       
       //mexPrintf("time0, sizes: %d, %d\n", grad_p1_xv_row, grad_p1_xv_col);

       //mexPrintf("1L\n");
       // this is grad_p1_v = ttt(grad_p1_x_v, f_mask(:,vi));
       // this can be further reduced using instead only the screened fasicles for vi
       // MARTHA: now we constantly update fasind (but not the considered atom indices), based on the level with which it decreases the objective
       // we only allow a small number of fascicles to become active for each voxel per iteration, and use a threshold improvement to stop these even sooner
       // fasind is updated within ttt_outer_product_spt, which updates the fasind for this voxel
       // assumes fasind has sufficient allocated space for the largest number of fascicles, so does not need reallocation
       begin = clock();
       if (flag_optimizefascicles == 1) {
           ttt_outer_product_spt(sub_grad_p1_v, val_grad_p1_v, &grad_p1_v_row, &grad_p1_v_col, sub_grad_p1_xv, val_grad_p1_xv, grad_p1_xv_row, grad_p1_xv_col, sub_fMask_v, val_fMask_v, fMask_v_row, fMask_v_col);
       }
       else {
           ttt_outer_product_spt(sub_grad_p1_v, val_grad_p1_v, &grad_p1_v_row, &grad_p1_v_col, sub_grad_p1_xv, val_grad_p1_xv, grad_p1_xv_row, grad_p1_xv_col, sub_fScreen_v, val_fScreen_v, fScreen_v_row, fScreen_v_col);
       }
       max_computed_size = max(grad_p1_v_row, max_computed_size);
       max_computed_sizes[1] = max(grad_p1_v_row, max_computed_sizes[1]);
       end = clock();
       time_spent[1] += (double)(end - begin) / CLOCKS_PER_SEC;
       
       //mexPrintf("time1, sizes: %d, %d\n", grad_p1_v_row, grad_p1_v_col);
       
       // the below is grad_p1_v = grad_p1_v .* entry_masks(:,vi,:);
       // this is not necessary since when calculating grad_p1_v, we already screen Phi
       // even more than entry_masks
       /**
       begin = clock();
       h_product_spt(sub_grad_p1_v_new, val_grad_p1_v_new, &grad_p1_v_new_row, &grad_p1_v_new_col, sub_grad_p1_v, val_grad_p1_v, grad_p1_v_row, grad_p1_v_col, sub_entryMask_v, val_entryMask_v, entryMask_v_row, entryMask_v_col, 0);
       max_computed_size = max(grad_p1_v_new_row, max_computed_size);
       max_computed_sizes[2] = max(grad_p1_v_new_row, max_computed_sizes[2]);
       end = clock();
       time_spent[2] += (double)(end - begin) / CLOCKS_PER_SEC;
       mexPrintf("time2, sizes: %d, %d\n", grad_p1_v_new_row, grad_p1_v_new_col);
       **/


       /************************ Finally sum up the three terms ******************/
       
       // this is g(:,vi,:) = g(:,vi,:) + grad_p1_v + lambda_g1 * grad_g1_v;
       
       sub_array = mxCreateCellArray(1, &sub_length);
       rc_pair_one = mxCreateDoubleMatrix(1, 1, mxREAL); 
       rc_pt_one = mxGetPr(rc_pair_one);
       rc_pair_two = mxCreateDoubleMatrix(1, 1, mxREAL); 
       rc_pt_two = mxGetPr(rc_pair_two);
       
       sub_for_g_v = mxCreateDoubleMatrix(default_eachG_row, default_eachG_col, mxREAL);
       sub_eachG = mxGetPr(sub_for_g_v);
       val_for_g_v = mxCreateDoubleMatrix(default_eachG_row, 1, mxREAL);
       val_eachG = mxGetPr(val_for_g_v);
       
       begin = clock();
       if (flag_optimizefascicles == 1) {
           add_subset(sub_eachG, val_eachG, &eachG_row, &eachG_col, sub_g_v, val_g_v, g_v_row, g_v_col, sub_grad_p1_v, val_grad_p1_v, grad_p1_v_row, grad_p1_v_col, sub_grad_g1_v, val_grad_g1_v, grad_g1_v_row, grad_g1_v_col, lambda_g1, sub_aScreen_v, sub_fMask_v, aScreen_v_row, fMask_v_row);
       }
       else {
           add_subset(sub_eachG, val_eachG, &eachG_row, &eachG_col, sub_g_v, val_g_v, g_v_row, g_v_col, sub_grad_p1_v, val_grad_p1_v, grad_p1_v_row, grad_p1_v_col, sub_grad_g1_v, val_grad_g1_v, grad_g1_v_row, grad_g1_v_col, lambda_g1, sub_aScreen_v, sub_fScreen_v, aScreen_v_row, fScreen_v_row);
       }
       max_computed_size = max(eachG_row, max_computed_size);
       max_computed_sizes[8] = max(eachG_row, max_computed_sizes[8]);
       end = clock();
       time_spent[8] += (double)(end - begin) / CLOCKS_PER_SEC;
       
       //mexPrintf("time8, sizes: %d, %d\n", eachG_row, eachG_col);
       //if (vi == 0)
           //print_array(sub_eachG, val_eachG, eachG_row, eachG_col);
       
       rc_pt_one[0] = eachG_row;
       rc_pt_two[0] = eachG_col;
       
       mxSetCell(sub_array, 0, sub_for_g_v);
       mxSetCell(sub_array, 1, val_for_g_v);
       mxSetCell(sub_array, 2, rc_pair_one);
       mxSetCell(sub_array, 3, rc_pair_two);
       mxSetCell(EachG_array, vi ,sub_array); 
        
       finalM += eachG_row;
   }
   if (max_computed_size > MAX_ROWS) {
       mexPrintf("Warning! MAX_ROWS %d must be bigger than any max_computed_size = %d\n", MAX_ROWS, max_computed_size);
   }
   mexPrintf("MAX_ROWS = %d must be bigger than any max_computed_size = %d\n", MAX_ROWS, max_computed_size);
   for (i = 0; i < 9; i++) {
       mexPrintf("%d = %d, time = %lf\n", i, max_computed_sizes[i], time_spent[i]);
   }

   plhs[0] = mxCreateDoubleMatrix(finalM,g_col, mxREAL); 
   double * sub = mxGetPr(plhs[0]); 
   plhs[1] = mxCreateDoubleMatrix(finalM, 1, mxREAL); 
   double * val = mxGetPr(plhs[1]); 
   
   int nonzero = 0;
   unsigned long m = 0;
   
   for(vi = 0; vi < Nvb; ++vi) {
       sub_array = mxGetCell(EachG_array, vi);
       sub_eachG = mxGetPr(mxGetCell(sub_array, 0));
       val_eachG = mxGetPr(mxGetCell(sub_array, 1));
       eachG_row = mxGetScalar(mxGetCell(sub_array, 2));
       eachG_col = mxGetScalar(mxGetCell(sub_array, 3));
      
       for(i = 0; i < eachG_row; ++i) {
           sub[get_index(m, 0, finalM, g_col)] = sub_eachG[get_index(i, 0, eachG_row, eachG_col)];
           sub[get_index(m, 1, finalM, g_col)] = vi + 1;
           sub[get_index(m, 2, finalM, g_col)] = sub_eachG[get_index(i, 1, eachG_row, eachG_col)];
           val[m] = val_eachG[i];
           ++m;
           if (sub_eachG[get_index(i, 0, eachG_row, eachG_col)] == 0 || sub_eachG[get_index(i, 1, eachG_row, eachG_col)] == 0) {
               nonzero++;
           }
       }
   }

   if (nonzero > 0) {
       mexPrintf("******************************Oh no! %d\n", nonzero);
   }
   mxFree(sub_grad_p1_xv);
   mxFree(val_grad_p1_xv);
   mxFree(sub_grad_p1_v);
   mxFree(val_grad_p1_v);
   mxFree(sub_grad_p1_v_new);
   mxFree(val_grad_p1_v_new);
   mxFree(sub_grad_g1_x);
   mxFree(val_grad_g1_x);
   mxFree(sub_grad_g1_x3);
   mxFree(val_grad_g1_x3);
   mxFree(sub_grad_g1_x4);
   mxFree(val_grad_g1_x4);
   mxFree(sub_temp);
   mxFree(val_temp);
   mxFree(sub_grad_g1_v);
   mxFree(val_grad_g1_v);
   //mxFree(sub_eachG);
   //mxFree(val_eachG);
   
   mxDestroyArray(aMask_array);
   mxDestroyArray(entryMask_array);   
   mxDestroyArray(fMask_array);
   mxDestroyArray(aScreen_array);
   mxDestroyArray(fScreen_array);
   mxDestroyArray(grad_p1_xt_array);
   mxDestroyArray(Phib_array);
   mxDestroyArray(g_array); 
   mxDestroyArray(Gv_array);
   mxDestroyArray(EachG_array);
   // subarray, sub_eachG are never initialized, they are just variables to refer to a part of EachG_array
} 