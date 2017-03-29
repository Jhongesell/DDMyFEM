//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema lineal usando factorizacion LU                            //
//                                                                                          //
// An�lisis y Dise�o y Programaci�n:                                                        //
//                                                                                          //
// Nombre:   Antonio Carrillo Ledesma                                                       //
// E-mail:   acl@www.mmc.igeofcu.unam.mx                                                    //
// P�gina:   http://www.mmc.igeofcu.unam.mx/acl                                             //
//                                                                                          //
//                                                                                          //
// Este programa es software libre. Puede redistribuirlo y/o modificarlo                    //
// bajo los t�rminos de la Licencia P�blica General de GNU seg�n es                         //
// publicada por la Free Software Foundation, bien de la versi�n 2 de                       //
// dicha Licencia o bien (seg�n su elecci�n) de cualquier versi�n                           //
// posterior.                                                                               //
//                                                                                          //
// Este programa se distribuye con la esperanza de que sea �til, pero SIN                   //
// NINGUNA GARANT�A, incluso sin la garant�a MERCANTIL impl�cita o sin                      //
// garantizar la CONVENIENCIA PARA UN PROP�SITO PARTICULAR. V�ase la                        //
// Licencia P�blica General de GNU para m�s detalles.                                       //
//                                                                                          //
// Deber�a haber recibido una copia de la Licencia P�blica General junto                    //
// con este programa. Si no ha sido as�, escriba a la Free Software                         //
// Foundation, Inc., en 675 Mass Ave, Cambridge, MA 02139, EEUU.                            //
//                                                                                          //
//                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////


#include "ResuelveFactorizacionLUDensa.hpp"


// Factoriza la matriz A en L y U dejandolas en la misma matriz
void ResuelveFactorizacionLUDensa::factoriza(void)
{
#ifdef DEPURAR
   if (!M->matrizCuadrada())
   {
      printf("\nError al factorizar por LU, matriz no cuadrada\n\n");
      exit(1);
   }
#endif

   int i, j, k;
   ldouble m, n;

   for (i = 0; i < M->renglones(); i++)
   {
#ifdef DEPURAR
      if (M->retorna(i,i) == 0.0)
      {
         printf("Error al factorizar el sistema en LU\n");
         exit(1);
      }
#endif
      for (j = i+1; j < M->renglones(); j++)
      {
         n = M->retorna(j,i);
         if (n == 0.0) continue;

         m = n / M->retorna(i,i);
         M->asigna(j, i, m);

         for (k = i+1; k < M->columnas(); k++)
         {
            M->asigna(j, k, M->retorna(j,k) - m * M->retorna(i,k) );
         }
      }
   }
   MatrizFactorizada = true;
}


// Resuelve el sistema lineal
void ResuelveFactorizacionLUDensa::resuelve(void)
{
   if (!MatrizFactorizada) factoriza();

   int i, j, n = M->renglones();
   ldouble t;

   Vector *y = new Vector(n);

   // Resuelve el sistema LY=B, dados L y B
   y->asigna(0,  B->retorna(0) / M->retorna(0,0) );
   for (i = 0; i < n; i++)
   {
      t = B->retorna(i);
      for (j = 0; j < i ; j++) t -= M->retorna(i,j) * y->retorna(j);
      y->asigna(i, t);
   }

   // Resuelve el sistema UX=Y, dados U  y Y
   X->asigna(n - 1, y->retorna(n - 1) / M->retorna(n - 1, n-1) );
   for (i = (n - 2); i >= 0; i--)
   {
      t = y->retorna(i);
      for (j = (n - 1); j > i; j--) t -= M->retorna(i,j) * X->retorna(j);
      X->asigna(i, t / M->retorna(i,j) );
   }

   delete y;
}
