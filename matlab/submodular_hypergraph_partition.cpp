#include <stdio.h>
#include <math.h>
#include <random>
#include <string>
#include <iostream>
#include "mex.h"
#include <Eigen/LU>

void swap_double(double * a, int p, int q)
{
    double temp = a[p];
    a[p] = a[q];
    a[q] = temp;
}

void copy(double * a, double * b, int len)
{
    // a <- b
    for(int i = 0; i < len; i++) { a[i] = b[i]; }
}

void all_ones_double(double * data, int len)
{
    for(int i = 0; i < len; i++) { data[i] = 1; }
}

double sum(double * data, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += data[i]; }
    return res;
}

double square_sum(double * data, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += data[i] * data[i]; }
    return res;
}

double inner(double * a, double * b, int len)
{
    double res = 0;
    for(int i = 0; i < len; i++) { res += a[i] * b[i]; }
    return res;
}

bool approx_eq(double a, double b)
{
    // check if a is close to b
    // ref: https://stackoverflow.com/questions/5595425
    double rtol = 1e-05;
    double atol = 1e-08;
    bool res = abs(a - b) <= fmax(rtol * fmax(abs(a), abs(b)), atol);
    return res;
}

bool approx_lte(double a, double b)
{
    // check if a <= b approximately
    double rtol = 1e-04; //1e-05;
    double atol = 1e-06; //1e-08;
    bool res = (a - b) <= fmax(rtol * fmax(abs(a), abs(b)), atol);
    return res;
}

void MYsort_descend(double * a, double * b, int len)
{
    // sort a, permutate b according to a
    if(len <= 1) return;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, len-1);
    int pos = distribution(generator), p = 0;
    double pivot = a[pos];
    swap_double(a, pos, len-1);
    swap_double(b, pos, len-1);
    for(int i = 0; i < len-1; i++)
    {
        if(a[i] > pivot)
        {
            swap_double(a, i, p);
            swap_double(b, i, p);
            p++;
        }
    }
    swap_double(a, p, len-1);
    swap_double(b, p, len-1);
    MYsort_descend(a, b, p);
    MYsort_descend(a+p+1, b+p+1, len-p-1);
}

void MYsort_ascend(double * a, double * b, int len)
{
    if(len <= 1) return;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, len-1);
    int pos = distribution(generator), p = 0;
    double pivot = a[pos];
    swap_double(a, pos, len-1);
    swap_double(b, pos, len-1);
    for(int i = 0; i < len-1; i++)
    {
        if(a[i] < pivot)
        {
            swap_double(a, i, p);
            swap_double(b, i, p);
            p++;
        }
    }
    swap_double(a, p, len-1);
    swap_double(b, p, len-1);
    MYsort_ascend(a, b, p);
    MYsort_ascend(a+p+1, b+p+1, len-p-1);
}

void remove_median(double * x, double * mu, int len)
{
    // return a non-constant vector x satisfying 0 ∈ argmin_c || x - c1 ||_{l1, mu}
    int i, pos=0;
    double * xsort = new double [len];
    double * mu_index = new double [len];
    copy(xsort, x, len);
    copy(mu_index, mu, len);
    MYsort_descend(xsort, mu_index, len);
    for(i = 1; i < len; i++)
        mu_index[i] += mu_index[i-1];
    for(i = 0; i < len; i++)
    {
        if(mu_index[i] >= mu_index[len-1]/2)
        {
            pos = i;
            break;
        }
    }
    for(i = 0; i < len; i++)
        x[i] -= xsort[pos];
    delete[] xsort;
    delete[] mu_index;
}

void edvw2para_e_part(double * subw, double * para, int len, int mode, double delta)
{
    // para has length len
    // subw has length len+1
    subw[0] = 0;    // the first entry
    subw[len] = 0;  // the last entry
    double total, ss;
    total = sum(para, len);
    for(int i = 1; i < len; i++)
    {
        ss = sum(para, i);
        switch(mode)
        {
            case 0: // c
                subw[i] = ss * (total - ss);
                break;
            case 1: // s
                subw[i] = std::min(ss, total - ss);
                break;
            case 2: // l
                subw[i] = std::min({ss, total - ss, delta * total});
                break;
            case 3: // ct
                subw[i] = std::min(ss * (total - ss), delta * total * total);
                break;
            default:
                throw std::invalid_argument( "specify a mode: {0: c, 1: s, 2: l, 3: ct}" );
        }
    }
}

void edmonds_greedy(double * q, double * x, double * para, int mode, double delta, double * a, int len)
{
    /*
      function
        - solve q ∈ argmin_{p ∈ B} <x, p>
        - B is the base polytope of some submodular function determined by para, mode, delta, a
      note
        - x, para, a have the same length len
    */
    double * index = new double [len];
    double * xsort = new double [len];
    double * psort = new double [len];
    double * asort = new double [len];
    double * subw = new double [len+1];
    for(int i = 0; i < len; i++)
    {
        index[i] = i;
        xsort[i] = x[i];
    }
    MYsort_ascend(xsort, index, len);
    for(int i = 0; i < len; i++)
    {
        psort[i] = para[(int)index[i]];
        asort[i] = a[(int)index[i]];
    }
    edvw2para_e_part(subw, psort, len, mode, delta);
    double temp = 0;
    for(int i = 0; i < len; i++)
    {
        temp += asort[i];
        subw[i+1] -= temp;
        q[(int)index[i]] = subw[i+1] - subw[i];
    }
    delete[] index;
    delete[] xsort;
    delete[] psort;
    delete[] asort;
    delete[] subw;
}

void affine_minimizer(double * arr_z, double * arr_a, const Eigen::MatrixXd &S, int nrow, int ncol)
{
    // solve min_{z ∈ aff(S)} || z ||
    // aff(S): the affine hall of the columns in S
    // return z and coefficients a
    if(ncol == 1) // S contains only one column
    {
        for(int i = 0; i < nrow; i++)
            arr_z[i] = S(i,0);
        arr_a[0] = 1;
    }
    else
    {
        Eigen::VectorXd z(nrow), a(ncol-1);
        Eigen::MatrixXd Q(nrow,ncol-1);
        Q = S.block(0,0,nrow,ncol-1).colwise() - S.col(ncol-1);
        a = (Q.transpose() * Q).inverse() * (-Q.transpose() * S.col(ncol-1)); // (ncol-1, 1)
        z = Q * a + S.col(ncol-1); // (nrow, 1)
        for(int i = 0; i < nrow; i++)
            arr_z[i] = z(i);
        for(int i = 0; i < ncol-1; i++)
            arr_a[i] = a(i);
        arr_a[ncol-1] = 1 - a.sum();
    }
}

void remove_col(Eigen::MatrixXd &matrix, int colToRemove, int numRows, int numCols)
{
    // ref: https://stackoverflow.com/a/46303314
    if(colToRemove < numCols)
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
    matrix.conservativeResize(numRows,numCols);
}

void remove_entry(double * data, int entryToRemove, int len)
{
    for(int i = entryToRemove; i < len-1; i++)
        data[i] = data[i+1];
}

void fujishige_wolfe_algorithm(double * output, double * a, double * para, int len, int mode, double delta, int maxiter_inner)
{
    /*
      find the point nearest to a in the base polytope of w_e(S)
      Fujishige-Wolfe algorithm:
        Compute the minimum norm point in the base polytope of some submodular function
    */

    double * start = new double [len];
    double * q = new double [len];
    double * x = new double [len];
    
    all_ones_double(start, len);
    edmonds_greedy(q, start, para, mode, delta, a, len);
    copy(x, q, len);
        
    int l = 1;
    int max_l = std::max(len, 1000);
    Eigen::MatrixXd S(len,1);
    for(int i = 0; i < len; i++) { S(i,0) = q[i]; } // set the first column to q
    double * lmbd = new double [max_l] {1}; // only consider the first l <= max_l entries
    
    bool flag;
    double theta, tmp;
    double * y = new double [len];
    double * coeff_extend = new double [max_l];
    
    // major cycle
    int cnt_major = 0;
    int cnt_minor = 0;
    while(1)
    {
        if(cnt_major > maxiter_inner){
            break;
        }
        cnt_major++;
        std::cout << "fw major: " << cnt_major << "; minor: " << cnt_minor << "\n";//"\r" << std::flush;
        // step (a)
        edmonds_greedy(q, x, para, mode, delta, a, len);
        // step (b)
        if(approx_lte(square_sum(x, len), inner(x, q, len)))
            break;
        // step (c)
        flag = false; // check if (updated) q has already existed in S
        for(int i = 0; i < l; i++) // for column
        {
            int cnt = 0;
            for(int j = 0; j < len; j++) // for row
            {
                if(approx_eq(q[j], S(j,i)))
                    cnt += 1;
            }
            if(cnt == len)
            {
                flag = true;
                break;
            }
        }
        if(!flag)
        {
            l += 1;
            S.conservativeResize(len, l);
            for(int i = 0; i < len; i++)
                S(i,l-1) = q[i]; // add a new column q
            lmbd[l-1] = 0;
        }
        // step (d) minor cycle
        cnt_minor = 0;
        while(1)
        {
            cnt_minor ++;
            // step (d-i)
            double * coeff = new double [l];
            affine_minimizer(y, coeff, S, len, l);
            // step (d-ii)
            flag = true; // check if coeff is all positive
            for(int i = 0; i < l; i++)
            {
                if(coeff[i] < 0)
                {
                    flag = false;
                    break;
                }
            }
            if(flag)
            {
                copy(coeff_extend, coeff, l);
                delete[] coeff;
                break;
            }
            // step (d-iii)
            theta = INFINITY;
            for(int i = 0; i < l; i++)
            {
                if(coeff[i] < 0)
                {
                    tmp = lmbd[i]/(lmbd[i] - coeff[i]);
                    theta = fmin(theta, tmp);
                }
            }

            for(int i = 0; i < len; i++)
                x[i] = theta * y[i] + (1 - theta) * x[i];
            for(int i = 0; i < l; i++)
                lmbd[i] = theta * coeff[i] + (1 - theta) * lmbd[i];
            for(int i = 0; i < l; i++)
            {
                if(approx_eq(lmbd[i], 0))
                {
                    remove_col(S, i, len, l-1);
                    remove_entry(lmbd, i, l);
                    l -= 1;
                }
            }
            copy(coeff_extend, coeff, l);
            delete[] coeff;
        }
        // step (e)
        copy(x, y, len);
        copy(lmbd, coeff_extend, l);
    }
    for(int i = 0; i < len; i++)
        output[i] = x[i] + a[i];
        
    delete[] start;
    delete[] q;
    delete[] x;
    delete[] lmbd;
    delete[] y;
    delete[] coeff_extend;
}

// int main()
// {
//     std::cout << "\n" << "j = " << "\n";
//     int len = 2;
//     double * a = new double [2]{-0.681667,-0.312701};
//     double * para = new double [2]{0.414049,0.445558};
//     double * output = new double [2];
//     for(int j = 0; j < 10000; j++)
//     {
//         std::cout << "\n" << "j = " << j << "\n";
//         fujishige_wolfe_algorithm(output, a, para, len, 0, 0);
//     }
//     return 0;
// }


double eval_Q1(double * x, int ** incidence_list, double ** parameter_list, int R, int * esize_list, int mode, double delta)
{
   double Q1 = 0;
   int esize;
   for(int i = 0; i < R; i++) // for each hyperedge
   {
       esize = esize_list[i];
       double * tempx = new double [esize];
       double * para = new double [esize];
       double * subw = new double [esize+1];
       for(int j = 0; j < esize; j++)
       {
           tempx[j] = x[incidence_list[i][j]];
           para[j] = parameter_list[i][j];
       }
       MYsort_descend(tempx, para, esize);
       edvw2para_e_part(subw, para, esize, mode, delta);
       for(int j = 0; j < esize; j++)
           Q1 += tempx[j] * (subw[j+1] - subw[j]);
       delete[] tempx;
       delete[] para;
       delete[] subw;
   }
   return Q1;
}

double eval_R1(double Q1, double * x, double * mu, int N)
{
   double temp = 0;
   for(int i = 0; i < N; i++)
       temp += abs(x[i]) * mu[i];
   double R1 = Q1/temp;
   return R1;
}

void derivative_mu_norm(double * g, double * x, double * mu, int N)
{
   double epsilon = 1e-10;
   double mu_x_pos = 0;
   double mu_x_neg = 0;
   double mu_x_zeros = 0;
   for(int i = 0; i < N; i++)
   {
       if(x[i] > epsilon)
           mu_x_pos += mu[i];
       else if(x[i] < -epsilon)
           mu_x_neg += mu[i];
       else
           mu_x_zeros += mu[i];
   }
   double muzero = (mu_x_neg - mu_x_pos)/mu_x_zeros;
   for(int i = 0; i < N; i++)
   {
       if(x[i] > epsilon)
           g[i] = mu[i];
       else if(x[i] < -epsilon)
           g[i] = -mu[i];
       else
           g[i] = muzero * mu[i];
   }
}

void optthreshold(double * labels, double * NCut,
                 double * x, int ** incidence_list, double ** parameter_list, double * mu,
                 int N, int R, int * esize_list, int mode, double delta)
{
   double * index = new double [N];
   double * xsort = new double [N];
   for(int i = 0; i < N; i++)
   {
       index[i] = i;
       xsort[i] = x[i];
   }
   MYsort_descend(xsort, index, N);
   int optpos = 0;
   double tempvol = 0, tempval, Q1;
   double * xtemp = new double [N](); // initialize all zeros
   *NCut = INFINITY;
   double sum_mu = sum(mu, N);
   for(int i = 0; i < N-1; i++)
   {
       xtemp[(int)index[i]] = 1;
       Q1 = eval_Q1(xtemp, incidence_list, parameter_list, R, esize_list, mode, delta);
       tempvol += mu[(int)index[i]];
       tempval = Q1/fmin(tempvol, sum_mu - tempvol);
       if(tempval < *NCut)
       {
           *NCut = tempval;
           optpos = i;
       }
   }
   for(int i = 0; i < N; i++)
   {
       if(i <= optpos)
           labels[(int)index[i]] = 1;
       else
           labels[(int)index[i]] = 0;
   }
   delete[] index;
   delete[] xsort;
   delete[] xtemp;
}

void main_func(double * labels, double * NCut, double * x_final,
       int ** incidence_list, double ** parameter_list, double * mu,
       int N, int R, int * esize_list, int mode, double delta,
       double dec_outloop, double err_inloop, double * warmstart,
       int maxiter_inner, int maxiter_outer)
{
   unsigned int T = 10000 * R;
   int record_dis = 100 * R;
   double Q1, tempeta, eta, normx, gap;
   int esize, picked;
   std::default_random_engine generator;
   std::uniform_int_distribution<int> distribution(0, R-1);

   double * x = new double [N];
   double ** a = new double * [R];
   double ** y = new double * [R];
   double ** newy = new double * [R];
   double * y_sum = new double [N]();
   double * g = new double [N];

   copy(x, warmstart, N);
   remove_median(x, mu, N);
   Q1 = eval_Q1(x, incidence_list, parameter_list, R, esize_list, mode, delta);
   tempeta = eval_R1(Q1, x, mu, N);

   for(int i = 0; i < R; i++)
   {
       esize = esize_list[i];
       a[i] = new double [esize]();
       y[i] = new double [esize];
       newy[i] = new double [esize];
       fujishige_wolfe_algorithm(y[i], a[i], parameter_list[i], esize, mode, delta, maxiter_inner);
       for(int j = 0; j < esize; j++)
           y_sum[incidence_list[i][j]] += y[i][j];
   }

   // step 1: repeat
   int cnt_outer = 0;
   while(1)
   {
       cnt_outer ++;
       std::cout << "\n a "<< cnt_outer << "\n";
       // step 2: update g
       derivative_mu_norm(g, x, mu, N);
       // step 3: RCDM
       eta = tempeta;
       for(int t = 0; t < T; t++)
       {
           picked = distribution(generator);
           for(int i = 0; i < esize_list[picked]; i++)
               a[picked][i] = y[picked][i] - y_sum[incidence_list[picked][i]] + eta * g[incidence_list[picked][i]];
           fujishige_wolfe_algorithm(newy[picked], a[picked], parameter_list[picked], esize_list[picked], mode, delta, maxiter_inner);
           for(int i = 0; i < esize_list[picked]; i++)
           {
               y_sum[incidence_list[picked][i]] += newy[picked][i] - y[picked][i];
               y[picked][i] = newy[picked][i];
           }
           if((t+1)%record_dis == 0)
           {
               std::cout << "-b";
               normx = 0;
               for(int i = 0; i < N; i++)
               {
                   x[i] = eta * g[i] - y_sum[i];
                   normx += x[i] * x[i];
               }
               normx = sqrt(normx);
               for(int i = 0; i < N; i++)
                   x[i] = x[i]/normx;
               Q1 = eval_Q1(x, incidence_list, parameter_list, R, esize_list, mode, delta);
               gap = Q1 - eta * inner(x, g, N) + normx;
               // step 4 & 5
               remove_median(x, mu, N);
               // step 6
               tempeta = eval_R1(Q1, x, mu, N);

               if(tempeta < (1-dec_outloop-1e-5)*eta || gap < err_inloop)
               {
                   if(tempeta < eta)
                       copy(x_final, x, N);
                   break;
               }
               
           }
       }
       if(tempeta > (1-dec_outloop)*eta)
           break;
       if(cnt_outer > maxiter_outer){
           break;
       }
   }

   optthreshold(labels, NCut, x_final, incidence_list, parameter_list, mu, N, R, esize_list, mode, delta);

    for (int i=0 ; i < R; i++){
        delete[] a[i];
        delete[] y[i];
        delete[] newy[i];
    }
    delete[] a;
    delete[] y;
    delete[] newy;
    delete[] x;
    delete[] y_sum;
    delete[] g;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs != 3 || nrhs != 12)
    {
        mexWarnMsgTxt("Check Parameters");
        return;
    }
    
    // INPUT
    
    double * mu = mxGetPr(prhs[2]);
    int N = *(mxGetPr(prhs[3]));
    int R = *(mxGetPr(prhs[4]));
    int mode = *(mxGetPr(prhs[5]));
    double delta = *(mxGetPr(prhs[6]));
    double dec_outloop = *(mxGetPr(prhs[7]));
    double err_inloop = *(mxGetPr(prhs[8]));
    double * warmstart = mxGetPr(prhs[9]);
    int maxiter_inner = *(mxGetPr(prhs[10]));
    int maxiter_outer = *(mxGetPr(prhs[11]));
    
    // read incidence_list
    const mxArray * incidence_list_org = prhs[0];
    mxArray * incidence_Element;
    int ** incidence_list = new int * [R];
    int * esize_list = new int [R];
    double * templist;
    for(int j = 0; j < R; j++)
    {
       incidence_Element = mxGetCell(incidence_list_org, j);
       esize_list[j] = (int)mxGetN(incidence_Element);
       incidence_list[j] = new int [esize_list[j]];
       templist = mxGetPr(incidence_Element);
       for(int k = 0; k < esize_list[j]; k++)
           incidence_list[j][k] = (int)templist[k]-1; // notice!
    }
    
    // read parameter_list
    const mxArray * parameter_list_org = prhs[1];
    mxArray * parameter_Element;
    double ** parameter_list = new double * [R];
    for(int j = 0; j < R; j++)
    {
       parameter_Element = mxGetCell(parameter_list_org, j);
       parameter_list[j] = new double [esize_list[j]];
       templist = mxGetPr(parameter_Element);
       for(int k = 0; k < esize_list[j]; k++)
           parameter_list[j][k] = templist[k];
    }

    // OUTPUT
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * labels = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double * NCut = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1, N, mxREAL);
    double * x_final = mxGetPr(plhs[2]);
    
    main_func(labels, NCut, x_final,
              incidence_list, parameter_list, mu,
              N, R, esize_list, mode, delta,
              dec_outloop, err_inloop, warmstart, 
              maxiter_inner, maxiter_outer);
   
}







