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



#ifndef __ResuelveJacobiDensa__
#define __ResuelveJacobiDensa__


#include "ResuelveJacobi.hpp"
#include "MatrizDensa.hpp"



/// Clase para resoluci�n del sistema lineal mediante Jacobi
/** @author Antonio Carrillo Ledesma
    @date primavera 2009
    @version 1.0.1
    @bug No hay errores conocidos
*/
class ResuelveJacobiDensa: public ResuelveJacobi
{

public:

   /// Constructor de la clase
   ResuelveJacobiDensa(void) : ResuelveJacobi()
   {
      M = NULL;
      RequiereMatriz = REQUIERE_MAT_DENS;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizDensa */
   ResuelveJacobiDensa(MatrizDensa *A) : ResuelveJacobi()
   {
      M = A;
      RequiereMatriz = REQUIERE_MAT_DENS;
   }

   /// Constructor de la clase
   /** @param A Puntero a una matriz del tipo MatrizDensa
       @param x Puntero a un Vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   ResuelveJacobiDensa(MatrizDensa *A, Vector *x, Vector *b) : ResuelveJacobi(x,b)
   {
      M = A;
      RequiereMatriz = REQUIERE_MAT_DENS;
   }

   /// Resuelve el sistema lineal
   void resuelve(void);

   /// Resuelve el sistema lineal
   /** @param x Puntero a un vector, soluci�n del sistema lineal
       @param b Puntero a un vector, lado derecho del sistema lineal */
   void resuelve(Vector *x, Vector *b)
   {
      X = x;
      B = b;
      resuelve();
   }

};

#endif
