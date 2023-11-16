    #include<stdio.h>

    #include<math.h>
    #include "InversionMatrice.h"
    
    
    double[N1+1] inv_mat(double a[N1+1,N1+1])

    {

      float k,d;
      
      k=N1+1

      int i, j;

      d = determinant(a, k);

      if (d == 0)

       printf("\nInverse of Entered Matrix is not possible\n");

      else

       cofactor(a, N1+1);

    }

     

    //Calcul le determinant de la matrice a

    float determinant(float a[N1+1][N1+1], float k)

    {

      float s = 1, det = 0, b[N1+1][N1+1];

      int i, j, m, n, c;

      if (k == 1)

        {

         return (a[0][0]);

        }

      else

        {

         det = 0;

         for (c = 0; c < k; c++)

           {

            m = 0;

            n = 0;

            for (i = 0;i < k; i++)

              {

                for (j = 0 ;j < k; j++)

                  {

                    b[i][j] = 0;

                    if (i != 0 && j != c)

                     {

                       b[m][n] = a[i][j];

                       if (n < (k - 2))

                        n++;

                       else

                        {

                         n = 0;

                         m++;

                         }

                       }

                   }

                 }

              det = det + s * (a[0][c] * determinant(b, k - 1));

              s = -1 * s;

              }

        }

     

        return (det);

    }

     
    
    void cofactor(float num[N1+1][N1+1], float f)

    {

     float b[N1+1][N1+1], fac[N1+1][N1+1];

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

                b[m][n] = num[i][j];

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

      transpose(num, fac, f);

    }

    //Trouve la transposé de la matrice

    void transpose(float num[N1+1][N1+1], float fac[N1+1][N1+1], float r)

    {

      int i, j;

      float b[N1+1][N1+1], inverse[N1+1][N1+1], d;

     

      for (i = 0;i < r; i++)

        {

         for (j = 0;j < r; j++)

           {

             b[i][j] = fac[j][i];

            }

        }

      d = determinant(num, r);

      for (i = 0;i < r; i++)

        {

         for (j = 0;j < r; j++)

           {

            inverse[i][j] = b[i][j] / d;

            }

        }
        

       return inverse

    }
