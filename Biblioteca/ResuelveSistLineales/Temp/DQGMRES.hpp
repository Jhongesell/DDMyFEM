//////////////////////////////////////////////////////////////////////////////////////////////
// Clase implementar DQGMRES con matrices bandadas o dispersas                              //
//                                                                                          //
// An�lisis y Dise�o y Programaci�n:                                                        //
//                                                                                          //
// Nombre:   Antonio Carrillo Ledesma                                                       //
// E-mail:   acl@www.mmc.geofisica.unam.mx                                                  //
// P�gina:   http://www.mmc.geofisica.unam.mx/acl                                           //
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



#ifndef __DQGMRES__
#define __DQGMRES__

#include <math.h>
#include <stdlib.h>
#include "MatrizBandDisp.hpp"

#define NMAXITER 50000

class DQGMRES 
{
protected:

   int n;             // vector size
   int k;             // # vectors to save
   int k1;
   int maxIter;       // maximum number of interations
   int nIter;         // current iteration
   //~ ArrayInt ind;         // Index array to p, Q
   //~ int index;
   double gm;         // gamma(m)
   double gm1;        // gamma(m + 1)
   double **p;      // [k + 1][n]
   double **cs;     // [k + 1][2] cosine, sine array for rotation
   double **h;      // Hessenberg
   double **q;      // quasi orthogonal vectors [k + 1][n]
   double *v;        // scratch
   double eps;
   int nMaxIter;

   /// Matriz Bandada o Dispersa
   MatrizBandDisp *M;


public:

   //~ DQGMRES(MultOp &mult, int k)
   //~ {
      //~ /* n = vector size, k = number of previous vectors to be saved & orthogonalized  */
      //~ this->eps = 1.0e-8;
      //~ this->nMaxIter = NMAXITER;
      //~ this->mult = &mult;
      //~ this->n = mult.getSize();
      //~ this->k = k;
      //~ inicializa();
   //~ }

   //~ DQGMRES(void)
   //~ {
      //~ this->eps = 1.0e-8;
      //~ this->nMaxIter = NMAXITER;
   //~ }

   /// Constructor de la clase
   DQGMRES(int n, MatrizBandDisp *M, int k) 
   {
      
      this->M = M;      
      this->n = n;
      this->k = k;
      this->eps = 1.0e-8;
      this->nMaxIter = NMAXITER;
      inicializa();
   }

   // Destructor de la clase
   ~DQGMRES()
   {
      delete M;
      libera();
   }
   
   /// vector size
   int getSize(void)
   {
      return n;
   }

   
   // Se quita pues genera un error en las referencias a MultOp
   //~ DQGMRES(MultOp &mult)
   //~ {
   //~ int n = (int) sqrt(mult.getSize());
   //~ DQGMRES(mult, n);
   //~ }

   void inicializa(void);

   void libera(void)
   {
      int i;
      for (i = 0; i < k1; i++)
      {
         delete []p[i];
         p[i] = NULL;
      }
      delete []p;
      p=NULL;

      for (i = 0; i < (k1+1); i++)
      {
         delete []q[i];
         q[i] = NULL;
      }
      delete []q;
      q=NULL;

      for (i = 0; i < k1; i++)
      {
         delete []cs[i];
         cs[i] = NULL;
      }
      delete []cs;
      cs=NULL;

      for (i = 0; i < (k1+1); i++)
      {
         delete []h[i];
         h[i] = NULL;
      }
      delete []h;
      h=NULL;

      delete []v;
   }

   void applyOmega(int m);

   void solve(Vector *x, Vector *b);

   int getIter(void)
   {
      return nIter;
   }

   void setMaxIter(int nmi)
   {
      nMaxIter = nmi;
   }

   void setEpsilon(double ep)
   {
      eps = ep;
   }

};

#endif
