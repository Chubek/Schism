#define SWAP(a, b) do { float __tmp = a; a = b; b = __tmp; } while(0)
#define SWAP_ROW(a, b, cnt)  do { for (int __i = 0; __i < cnt; __i++)  { SWAP(a[__i], b[__i]); } } while(0)
#define COPY_VEC(a, b, cnt)  do { for (int __iii = 0; __iii < cnt; __iii++) { b[__iii] = a[__iii]; } } while(0)
#define ELEM_WISE_VEC_VEC(a, b, cnt) do { for (int __ii = 0; __ii < cnt; __ii++) { a[__ii] -= b[__ii]; } } while(0)
#define ELEM_WISE_VEC_SCALAR(a, b, cnt) do { for (int __j = 0; __j < cnt; __j++ ) { a[__j] *= b; } } while(0)

// P = L = I
// U = A; main matrix
// A = [count][count] square matrix
extern "C" __global__ void lu_decompose_pivot(float** L, float** U, float** P, int count) {
    for (int t = blockIdx.x * blockDim.x + threadIdx.x; t < count; t += blockDim.x * gridDim.x) {
        for (int i = 0; i < t; i++) {
            int k = i;

            while (U[i][i] == 0) {
                SWAP_ROW(U[i], U[k + 1], count);
                SWAP_ROW(P[i], P[k + 1], count);

                ++k;
            }

            for (int j = i + 1; j < t; j++) {
                L[j][i] = U[j][i] / U[i][i];

                float cpy[count];
                COPY_VEC(U[i], cpy, count);

                ELEM_WISE_VEC_SCALAR(cpy, L[j][i], count);
                ELEM_WISE_VEC_VEC(U[j], cpy, count);
            }


        }
    }
}