//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para resolver un sistema lineal mediante Jacobi                                    //
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


#include <math.h>
#include "ResuelveJacobiDensa.hpp"


// Resuelve el sistema lineal
void ResuelveJacobiDensa::resuelve(void)
{
#ifdef DEPURAR
   if (!M->matrizCuadrada()) error("matriz no cuadrada");
#endif

   int i, j, sw;
   ldouble sum;
   Vector *xt = new Vector(M->renglones());
   if (!xt) error("no hay memoria para el vector auxiliar");

   for (NumIteraciones = 1; NumIteraciones <= Iter; NumIteraciones ++)
   {
      for (i = 0; i < M->renglones(); i++)
      {
         if (M->retorna(i,i) == 0.0)
         {
            delete xt;
            return;
         }
         sum = 0.0;
         for (j = 0; j < M->renglones(); j ++)
         {
            if (i == j) continue;
            sum += M->retorna(i,j) * xt->retorna(j);
         }
         xt->asigna(i,(B->retorna(i) - sum) / M->retorna(i,i));
      }
      sw = 0;
      for (i = 0; i < M->renglones(); i++)
      {
#ifdef __Double__
         if (fabs(X->retorna(i) - xt->retorna(i)) > Ep) sw = 1;
#else
         if (fabsl(X->retorna(i) - xt->retorna(i)) > Ep) sw = 1;
#endif
         X->asigna(i,xt->retorna(i));
      }
      if (!sw) break;
   }
   delete xt;
}
