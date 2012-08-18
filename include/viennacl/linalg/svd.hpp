#ifndef VIENNACL_LINALG_SVD_HPP
#define VIENNACL_LINALG_SVD_HPP

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/svd.hpp
    @brief Provides singular value decomposition using a block-based approach.  Experimental in 1.3.x.
    
    Contributed by Volodymyr Kysenko.
*/


// Note: Boost.uBLAS is required at the moment
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


#include <cmath>

#include "viennacl/matrix.hpp"
#include "viennacl/linalg/kernels/svd_kernels.h"

namespace viennacl 
{
  namespace linalg 
  {
  
    //const std::string SVD_KERNELS_FOLDER = "../../non-release/svd-kernels/";
    //const std::string SVD_BIDIAG_PROGRAM = "bidiag.cl";

    const std::string SVD_BIDIAG_PACK_KERNEL = "bidiag_pack";
    const std::string SVD_HOUSEHOLDER_COL_KERNEL = "house_col";
    const std::string SVD_HOUSEHOLDER_ROW_KERNEL = "house_row";
    const std::string SVD_COPY_COL_KERNEL = "copy_col";
    const std::string SVD_COPY_ROW_KERNEL = "copy_row";
    const std::string SVD_MATRIX_TRANSPOSE_KERNEL = "transpose_inplace";
    const std::string SVD_INVERSE_SIGNS_KERNEL = "inverse_signs";
    const std::string SVD_GIVENS_PREV_KERNEL = "givens_prev";
    
    namespace detail 
    {
      static const float EPS = 0.00001f;
      static const std::size_t ITER_MAX = 50;

      inline float pythag(float a, float b) 
      {
        float absa = std::abs(a);
        float absb = std::abs(b);

        if(absa > absb) {
          return absa * sqrt(1.0f + pow(absb / absa, 2));
        } else {
          return absb * sqrt(1.0f + pow(absa / absb, 2));
        }
      }

      inline float sign(float val) 
      {
          return val >= 0.0f ? 1.0f : -1.0f;
      }

      inline float norm_lcl(std::vector<float>& x, unsigned int size) 
      {
        float x_norm = 0.0;
        for(std::size_t i = 0; i < size; i++) x_norm += std::pow(x[i], 2);
        x_norm = std::sqrt(x_norm);
        return x_norm;
      }

      template <typename T>
      void normalize(std::vector<T>& x, unsigned int size) 
      {
        float x_norm = norm_lcl(x, size);
        for(std::size_t i = 0; i < size; i++) {
            x[i] /= x_norm;
        }
      }

      template <typename T>
      void householder_vector(std::vector<T> & v, unsigned int start)
      {
        float x_norm = norm_lcl(v, v.size());
        float alpha = -sign(v[start]) * x_norm;
        v[start] += alpha;
        normalize(v, v.size());
      }

      template <typename MatrixType>
      void transpose(MatrixType& A)
      {

        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(), SVD_MATRIX_TRANSPOSE_KERNEL);

        viennacl::ocl::enqueue(kernel(
                                      A,
                                      static_cast<cl_uint>(A.internal_size1()),
                                      static_cast<cl_uint>(A.internal_size2())
                              ));
      }

      template<typename MatrixType, typename VectorType>
      void givens_prev(MatrixType& matrix,
                        VectorType& tmp1,
                        VectorType& tmp2,
                        int n,
                        int l,
                        int k
                      )
      {
        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(), SVD_GIVENS_PREV_KERNEL);

        kernel.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(viennacl::traits::size1(matrix), 256));
        kernel.local_work_size(0, 256);

        viennacl::ocl::enqueue(kernel(
                                      matrix,
                                      tmp1,
                                      tmp2,
                                      static_cast<cl_uint>(n),
                                      static_cast<cl_uint>(matrix.internal_size1()),
                                      static_cast<cl_uint>(l + 1),
                                      static_cast<cl_uint>(k + 1)
                              ));
      }


      template<typename MatrixType, typename VectorType>
      void change_signs(MatrixType& matrix, VectorType& signs, int n)
      {
        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(), SVD_INVERSE_SIGNS_KERNEL);

        kernel.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(viennacl::traits::size1(matrix), 16));
        kernel.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(viennacl::traits::size2(matrix), 16));

        kernel.local_work_size(0, 16);
        kernel.local_work_size(1, 16);

        viennacl::ocl::enqueue(kernel(
                                      matrix,
                                      signs,
                                      static_cast<cl_uint>(n),
                                      static_cast<cl_uint>(matrix.internal_size1())
                              ));
      }

      template<typename MatrixType>
      void svd_qr_shift(MatrixType& vcl_u,
                        MatrixType& vcl_v,
                        boost::numeric::ublas::vector<float> &q, 
                        boost::numeric::ublas::vector<float> &e)
      {
        int n = q.size();
        int m = vcl_u.size1();

        detail::transpose(vcl_u);
        detail::transpose(vcl_v);

        std::vector<float> signs_v(n, 1.0f);
        std::vector<float> cs1(n), ss1(n), cs2(n), ss2(n);
        
        viennacl::vector<float> tmp1(n), tmp2(n);

        bool goto_test_conv = false;

        for (int k = n - 1; k >= 0; k--) {
          // std::cout << "K = " << k << std::endl;

          std::size_t iter = 0;
          for (iter = 0; iter < detail::ITER_MAX; iter++) {
            // test for split
            int l;
            for (l = k; l >= 0; l--) {
              goto_test_conv = false;
              if (fabs(e[l]) <= detail::EPS) {
                // set it
                goto_test_conv = true;
                break;
              }

              if (fabs(q[l - 1]) <= detail::EPS) {
                // goto
                break;
              }
            }

            if (!goto_test_conv) {
              float c = 0.0;
              float s = 1.0;

              //int l1 = l - 1;
              int l2 = k;

              for (int i = l; i <= k; i++) {
                float f = s * e[i];
                e[i] = c * e[i];

                if (fabs(f) <= detail::EPS) {
                  l2 = i - 1;
                  break;
                }

                float g = q[i];
                float h = detail::pythag(f, g);
                q[i] = h;
                c = g / h;
                s = -f / h;

                cs1[i] = c;
                ss1[i] = s;
              }

              // std::cout << "Hitted!" << l1 << " " << l2 << "\n";

              // for(int i = l; i <= l2; i++) 
              // {
              //   for (int j = 0; j < m; j++) 
              //   {
              //     float y = u(j, l1);
              //     float z = u(j, i);
              //     u(j, l1) = y * cs1[i] + z * ss1[i];
              //     u(j, i) = -y * ss1[i] + z * cs1[i];
              //   }
              // }
            }

            float z = q[k];

            if (l == k) {
              if (z < 0.0f) {
                q[k] = -z;

                signs_v[k] *= -1.0f;
              }

              break;
            }

            if (iter >= detail::ITER_MAX - 1) {
              break;
            }

            float x = q[l];
            float y = q[k - 1];
            float g = e[k - 1];
            float h = e[k];
            float f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
            
            g = detail::pythag(f, 1.0);

            if (f < 0) {
              f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
            } else {
              f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
            }

            float c = 1.0;
            float s = 1.0;

            for (std::size_t i = l + 1; i <= static_cast<std::size_t>(k); i++) 
            {
              g = e[i];
              y = q[i];
              h = s * g;
              g = c * g;
              float z = detail::pythag(f, h);
              e[i - 1] = z;
              c = f / z;
              s = h / z;
              f = x * c + g * s;
              g = -x * s + g * c;
              h = y * s;
              y = y * c;
              
              cs1[i] = c;
              ss1[i] = s;

              z = detail::pythag(f, h);
              q[i - 1] = z;
              c = f / z;
              s = h / z;
              f = c * g + s * y;
              x = -s * g + c * y;

              cs2[i] = c;
              ss2[i] = s;
            }
            
            {
              viennacl::copy(cs1, tmp1);
              viennacl::copy(ss1, tmp2);

              givens_prev(vcl_v, tmp1, tmp2, n, l, k);
            }

            {
              viennacl::copy(cs2, tmp1);
              viennacl::copy(ss2, tmp2);

              givens_prev(vcl_u, tmp1, tmp2, m, l, k);
            }
            
            e[l] = 0.0;
            e[k] = f;
            q[k] = x;
          }

        }

        
        viennacl::copy(signs_v, tmp1);
        change_signs(vcl_v, tmp1, n);

        // transpose singular matrices again
        detail::transpose(vcl_u);
        detail::transpose(vcl_v);
      }

      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      void eye(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A)
      {
      
        std::vector<SCALARTYPE> foo(A.size1() * A.size1(), 0);
        
        for(std::size_t i = 0; i < A.size1(); i++)
        {
          foo[i*A.size1() + i] = 1;
        }

        viennacl::fast_copy(&foo[0], &foo[0] + foo.size(), A);
      }
      
      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      void copy_vec(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                    viennacl::vector<SCALARTYPE, ALIGNMENT>& V,
                    std::size_t row_start, 
                    std::size_t col_start, 
                    bool copy_col
      )
      {

        std::string kernel_name = copy_col ? SVD_COPY_COL_KERNEL : SVD_COPY_ROW_KERNEL;
        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(),
                                                                  kernel_name);

        viennacl::ocl::enqueue(kernel(
                                      A, 
                                      V, 
                                      static_cast<cl_uint>(row_start), 
                                      static_cast<cl_uint>(col_start),
                                      copy_col ? static_cast<cl_uint>(A.size1())
                                               : static_cast<cl_uint>(A.size2()),
                                      static_cast<cl_uint>(A.internal_size2())
                              ));

      }

      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      bool householder_c(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                          viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& Q,
                          viennacl::vector<SCALARTYPE, ALIGNMENT>& D,
                          std::size_t start) 
      {

        std::size_t row_start = start, col_start = start;

        if(row_start + 1 >= A.size1()) 
          return false;

        std::vector<float> tmp(A.size1(), 0);

        copy_vec(A, D, row_start, col_start, true);
        fast_copy(D.begin(), D.begin() + (A.size1() - row_start), tmp.begin() + row_start);

        detail::householder_vector(tmp, row_start);
        fast_copy(tmp, D);

        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(), SVD_HOUSEHOLDER_COL_KERNEL);

        //kernel.global_work_size(0, A.size1() << 1);

        viennacl::ocl::enqueue(kernel(
                                      A,
                                      Q,
                                      D,
                                      static_cast<cl_uint>(row_start),
                                      static_cast<cl_uint>(col_start),
                                      static_cast<cl_uint>(A.size1()),
                                      static_cast<cl_uint>(A.size2()),
                                      static_cast<cl_uint>(A.internal_size2()),
                                      static_cast<cl_uint>(Q.internal_size2()),
                                      viennacl::ocl::local_mem(static_cast<cl_uint>(128 * 4))
                              ));

        return true;
      }

      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      bool householder_r(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                          viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& Q,
                          viennacl::vector<SCALARTYPE, ALIGNMENT>& S,
                          std::size_t start)
      {
      
        std::size_t row_start = start, col_start = start + 1;
        if(col_start + 1 >= A.size2()) 
          return false;

        std::vector<float> tmp(A.size2(), 0);

        copy_vec(A, S, row_start, col_start, false);
        fast_copy(S.begin(), S.begin() + (A.size2() - col_start), tmp.begin() + col_start);

        detail::householder_vector(tmp, col_start);
        fast_copy(tmp, S);

        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(), SVD_HOUSEHOLDER_ROW_KERNEL);

        viennacl::ocl::enqueue(kernel(
                                      A,
                                      Q,
                                      S,
                                      static_cast<cl_uint>(row_start),
                                      static_cast<cl_uint>(col_start),
                                      static_cast<cl_uint>(A.size1()),
                                      static_cast<cl_uint>(A.size2()),
                                      static_cast<cl_uint>(A.internal_size2()),
                                      static_cast<cl_uint>(Q.internal_size2()),
                                      viennacl::ocl::local_mem(static_cast<cl_uint>(128 * 4))
                                ));
        return true;
      }

      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      void bidiag_pack(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                        viennacl::vector<SCALARTYPE, ALIGNMENT>& D,
                        viennacl::vector<SCALARTYPE, ALIGNMENT>& S
                      )
      {
        viennacl::ocl::kernel& kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::svd<float, 1>::program_name(), SVD_BIDIAG_PACK_KERNEL);

        viennacl::ocl::enqueue(kernel(
                                      A, 
                                      D, 
                                      S,
                                      static_cast<cl_uint>(A.size1()), 
                                      static_cast<cl_uint>(A.size2()),
                                      static_cast<cl_uint>(A.internal_size2())
                                    ));
      }

      template <typename SCALARTYPE, unsigned int ALIGNMENT>
      void bidiag(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& Ai,
                  viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& QL,
                  viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& QR)
      {
        std::size_t row_num = Ai.size1();
        std::size_t col_num = Ai.size2();

        std::size_t to = std::min(row_num, col_num);
        std::size_t big_to = std::max(row_num, col_num);

        //for storing householder vector
        viennacl::vector<SCALARTYPE, ALIGNMENT> hh_vector(big_to);

        eye(QL);
        eye(QR);

        for(std::size_t i = 0; i < to; i++) 
        {
          householder_c(Ai, QL, hh_vector, i);
          householder_r(Ai, QR, hh_vector, i);
        }
      }

    } // namespace detail


    /** @brief Computes the singular value decomposition of a matrix A. Experimental - works for single precision (float) only. Experimental in 1.3.x
     * 
     * @param A     The input matrix. Will be overwritten with a diagonal matrix containing the singular values on return
     * @param QL    The left orthogonal matrix
     * @param QR    The right orthogonal matrix
     */
    template <unsigned int ALIGNMENT>
    void svd(viennacl::matrix<float, row_major, ALIGNMENT> & A,
              viennacl::matrix<float, row_major, ALIGNMENT> & QL,
              viennacl::matrix<float, row_major, ALIGNMENT> & QR) 
    {
      typedef float SCALARTYPE;
      
      viennacl::linalg::kernels::svd<SCALARTYPE, 1>::init();

      std::size_t row_num = A.size1();
      std::size_t col_num = A.size2();

      std::size_t to = std::min(row_num, col_num);


      viennacl::vector<SCALARTYPE, ALIGNMENT> d(to);
      viennacl::vector<SCALARTYPE, ALIGNMENT> s(to + 1);
      
      // first stage
      detail::bidiag(A, QL, QR);
      detail::bidiag_pack(A, d, s);

      // second stage
      boost::numeric::ublas::vector<SCALARTYPE> dh(to, 0.0f);
      boost::numeric::ublas::vector<SCALARTYPE> sh(to + 1, 0.0f);

      boost::numeric::ublas::matrix<float> h_U(row_num, row_num);
      boost::numeric::ublas::matrix<float> h_V(col_num, col_num);

      fast_copy(d, dh);
      fast_copy(s, sh);

      detail::svd_qr_shift( QL, QR, dh, sh);

      boost::numeric::ublas::matrix<float> h_Sigma(row_num, col_num);
      h_Sigma.clear();

      for (std::size_t i = 0; i < to; i++)
        h_Sigma(i, i) = dh(i);

      copy(h_Sigma, A);
    }
  }
}
#endif
