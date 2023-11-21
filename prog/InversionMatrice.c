    #include<stdio.h>

    #include<math.h>
    #include "InversionMatrice.h"

    
    void inv_mat(double a[N1+1][N1+1])

    {

      double d;
      
      int i, j;

      d = determinant(a, N1+1);

      if (d == 0)

       printf("\nInverse of Entered Matrix is not possible\n");

      else

       cofactor(a, N1+1);

    }

     

    /*For calculating Determinant of the Matrix */

    double determinant(double a[N1+1][N1+1], double k1)

    {

      double s = 1, det = 0, b[N1+1][N1+1];

      int i, j, m, n, c;

      if (k1 == 1)

        {

         return (a[0][0]);

        }

      else

        {

         det = 0;

         for (c = 0; c < k1; c++)

           {

            m = 0;

            n = 0;

            for (i = 0;i < k1; i++)

              {

                for (j = 0 ;j < k1; j++)

                  {

                    b[i][j] = 0;

                    if (i != 0 && j != c)

                     {

                       b[m][n] = a[i][j];

                       if (n < (k1 - 2))

                        n++;

                       else

                        {

                         n = 0;

                         m++;

                         }

                       }

                   }

                 }

              det = det + s * (a[0][c] * determinant(b, k1 - 1));

              s = -1 * s;

              }

        }

     

        return (det);

    }

     

    void cofactor(double num1[N1+1][N1+1], double f)

    {

     double b[N1+1][N1+1], fac[N1+1][N1+1];

     int p, q, m, n, i, j;

     for (q = 0;q < f; q++)

     {

       for (p = 0;p < f; p++)

        {

         m = 0;

         n = 0;

         for (i = 0;i < f; i++)

         {

           for (j = 0;j < f; j++)

            {

              if (i != q && j != p)

              {

                b[m][n] = num1[i][j];

                if (n < (f - 2))

                 n++;

                else

                 {

                   n = 0;

                   m++;

                   }

                }

            }

          }

          fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);

        }

      }

      transpose(num1, fac, f);

    }

    /*Finding transpose of matrix*/ 

    void transpose(double num1[N1+1][N1+1], double fac[N1+1][N1+1], double r)

    {

      int i, j;

      double b[N1+1][N1+1], inverse[N1+1][N1+1], d;

     

      for (i = 0;i < r; i++)

        {

         for (j = 0;j < r; j++)

           {

             b[i][j] = fac[j][i];

            }

        }

      d = determinant(num1, r);

      for (i = 0;i < r; i++)

        {

         for (j = 0;j < r; j++)

           {

            inverse[i][j] = b[i][j] / d;

            }

        }
    

    }
