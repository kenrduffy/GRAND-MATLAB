/* Basic ORBGRAND (with list output)
 *
 * Install: mex -O orbGRAND.c
 *
 * [chat, LW, NT] = orbGRAND(llr, H, L, Tmax);
 *
 * n     :                 codelength
 *
 * chat  :(n,L) matrix,    L best candidates c
 * score :(1,L) matrix,    Logistic weight of the candidates
 * NT    :scalar,          Number of guesses
 *
 * LLR   :(n,1) matrix,    log likelihood ratio of the channel outputs
 * H     :(n*s,1) matrix,    Parity-check matrix
 * L     :scalar,          list size of the list decoder
 * Tmax  :scalar,          Maximum number of guesses
 *
 * References:
 * [1] K. R. Duffy, W. An, and M. Médard. "Ordered reliability bits guessing random additive noise decoding." IEEE Transactions on Signal Processing 70 (2022): 4528-4542.
 *
 * 13 Jan 2023
 */

/* C Functions */
#include <mex.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#define Inf 0x7fffffff

/* Parity_Check */
uint8_t ParityCheck(uint8_t *c, uint8_t *H, uint64_t n, uint64_t s) {
    uint8_t syndrome;
    for (size_t j = 0; j < s; j++){
        syndrome = 0;
        for (size_t i = 0; i < n; i++)
            syndrome ^= ( c[i] * H[j*n + i] );
        if (syndrome == 1)
            return 0;
    }
    return 1;
}

/* Hard Decision */
void HardDec(uint8_t *c, double *llr, uint64_t n) {
    for (size_t i = 0; i < n; i++){
        if (llr[i] > 0.0)
            c[i] = 0;
        else
            c[i] = 1;
    }
}

/* Add TEP */
void AddTEP(uint8_t *c, uint8_t *cHD, uint8_t *TEP, size_t *perm, uint64_t n) {
    for (size_t i = 0; i < n; i++)
        c[perm[i]] = cHD[perm[i]] ^ TEP[i];
}

/* Quick Sort (ascend) */
void QuickSort (double *a, size_t *perm, uint64_t n) {
    if (n < 2) return;
    uint64_t i, j;
    double t, p;
    size_t tt;
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p) i++;
        while (p < a[j]) j--;
        if (i >= j) break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
        tt = perm[i];
        perm[i] = perm[j];
        perm[j] = tt;
    }
    QuickSort(a, perm, i);
    QuickSort(a + i, perm + i, n - i);
}

/* Logistic Weight */
uint64_t findLW(uint8_t *TEP, size_t a, size_t b) {
    uint64_t LW = 0;
    for (size_t i = a; i < b+1; i++){
        if (TEP[i-1] == 1)
            LW += i;
    }
    return LW;
}

uint8_t isLastTEP(uint8_t *TEP, uint64_t n) {
    uint64_t LW = findLW(TEP, 1, n);
    if (LW == 0)
        return 1;
    for (size_t i = n; i > 0; i--){
        if ( (i-1)*i < 2*LW ){
            if (TEP[i-1] == 0)
                return 0;
            LW -= i;
            if (LW == 0)
                return 1;
        }
    }
    return 0;
}

void findLastTEP(uint8_t *TEP, uint64_t LW, uint64_t n) {
    for (size_t i = 0; i < n; i++)
        TEP[i] = 0;
    for (size_t i = n; i > 0; i--){
        if ( i <= LW ){
            TEP[i-1] = 1;
            LW -= i;
            if (LW == 0)
                return;
        }
    }
}

void NextTEP(uint8_t *TEP, uint64_t n) {
    size_t LeftMost = 0;
    uint64_t r;
    uint64_t LW = findLW(TEP, 1, n);
    if (isLastTEP(TEP, n) == 1){
        findLastTEP(TEP, LW+1, n);
        return;
    }else{
        while(true){
            for (size_t i = 3; i < n+1; i++){
                if (TEP[i-1] == 1){
                    LeftMost = i;
                    break;
                }
            }
            r = findLW(TEP, LeftMost+1, n);
            if ( ((LeftMost-1)*LeftMost + 2*r) < 2*LW ){
                for (size_t i = 1; i < LeftMost+1; i++)
                    TEP[i-1] = 0;
            }else{
                findLastTEP(TEP, LW-r, LeftMost-1);
                TEP[LeftMost-1] = 0;
                return;
            }
        }
    }
    return;
}

/* Main Function */
void orbGRAND(double *chat, double *score, double *T, double *llr, uint8_t *H, uint64_t n, uint64_t s, uint64_t L, uint64_t Tmax){
    /* Create vectors */
    size_t *perm = calloc(n, sizeof(size_t));
    for(size_t i = 0; i < n; i++)
        perm[i] = i;
    uint8_t *cHD = calloc(n, sizeof(uint8_t));
    uint8_t *TEP = calloc(n, sizeof(uint8_t));
    uint8_t *c   = calloc(n, sizeof(uint8_t));
    double *absL = calloc(n, sizeof(double));
    for(size_t i = 0; i < L; i++)
        score[i] = Inf;
    uint64_t cur_L = 0;
    /* Initialize */
    HardDec(cHD, llr, n);
    T[0] = 1;
    
    if (ParityCheck(cHD, H, n, s) == 1){
        score[cur_L] = 0;
        for(size_t i = 0; i < n; i++)
            chat[cur_L*n + i] = cHD[i];
        cur_L++;
    }
    /* Main GRAND */
    for (size_t i = 0; i < n; i++){
        TEP[i] = 0;
        absL[i] = fabs(llr[i]);
    }
    QuickSort (absL, perm, n);
    while ( (cur_L < L) && (T[0] < Tmax) ){
        NextTEP(TEP, n);
        AddTEP(c, cHD, TEP, perm, n);
        T[0]++;
        if (ParityCheck(c, H, n, s) == 1){
            score[cur_L] = findLW(TEP, 1, n);
            for(size_t i = 0; i < n; i++)
                chat[cur_L*n + i] = c[i];
            cur_L++;
        }
    }
    
    
    /* Clean up allocated memory */
    free(perm);
    free(cHD);
    free(TEP);
    free(c);
    free(absL);
}

/* Mexfunction Interface */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    /* Check for proper number of arguments */
    if( nrhs != 4 )
        mexErrMsgTxt("Must have 3 input arguments.");
    if( nlhs != 3 )
        mexErrMsgTxt("Must have 3 output arguments.");
    if( mxGetN(prhs[0]) != 1 || !mxIsClass(prhs[0], "double"))
        mexErrMsgTxt("First Input (LLR) must be a column-vector of type double.");
    uint64_t n = mxGetNumberOfElements(prhs[0]); // code length n
    if( mxGetN(prhs[1]) != 1 || !mxIsClass(prhs[1], "uint8"))
        mexErrMsgTxt("Second Input (Parity check matrix) must be a column-vector of length s*n and type uint8.");
    uint64_t N = mxGetNumberOfElements(prhs[1]); // s*n
    uint64_t s = (uint64_t) N/n; // number of constraints
    if( mxGetNumberOfElements(prhs[2]) != 1 || !mxIsClass(prhs[2], "uint64"))
        mexErrMsgTxt("Third Input (List size) must be an uint64 scalar.");
    if( mxGetNumberOfElements(prhs[3]) != 1 || !mxIsClass(prhs[3], "uint64"))
        mexErrMsgTxt("Forth Input (Maximum number of guesses) must be an uint64 scalar.");
    /* input */
    double *llr = mxGetPr(prhs[0]);
    uint8_t *H = (uint8_t *)mxGetData(prhs[1]);
    uint64_t L = (uint64_t)mxGetScalar(prhs[2]);
    uint64_t Tmax = (uint64_t)mxGetScalar(prhs[3]);
    /* create the output vector */
    plhs[0] = mxCreateDoubleMatrix(n, L, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, L, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *chat  = mxGetData(plhs[0]);
    double *score = mxGetData(plhs[1]);
    double *T     = mxGetData(plhs[2]);
    /* use C functions in mexfunction */
    orbGRAND(chat, score, T, llr, H, n, s, L, Tmax);
}