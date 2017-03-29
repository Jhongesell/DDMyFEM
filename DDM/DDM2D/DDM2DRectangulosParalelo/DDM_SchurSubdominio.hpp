//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para definir los comportamientos de un subdominio mediante el m�todo de descompo-  //
// sici�n de dominio de subestructuraci�n Schur                                             //
//                                                                                          //
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



#ifndef __DDM_SchurSubdominio__
#define __DDM_SchurSubdominio__



#include "VectorExtendido.hpp"
#include "MatrizBandDisp.hpp"
#include "ResuelveSistemaLineal.hpp"
#include "Geometria.hpp"
#include "CtrlMemoria.hpp"
#include "Problema_2D.hpp"
#include "Geometria_2DRectangulos.hpp"




// Usar Factorizacion LU o CGM para resolver la proyeccion de la frontera interior
//~ #define USANDO_FACTORIZACION_LU




/// Clase para definir los comportamientos de un subdominio mediante el m�todo de descomposici�n de dominio de subestructuraci�n (Complemento de Schur)
/** @author Antonio Carrillo Ledesma
    @date verano 2009
    @version 1.0.0
    @bug No hay errores conocidos
*/
class DDM_SchurSubdominio
{

protected:

   /// Puntero a la definicion de la geometria
   Geometria      *pGe;
   /// Puntero a la definicion del problema
   Problema_2D    *pPr2D;
   /// Puntero a la definicion de la geometria
   Geometria_2D   *pGe2D;
   /// Objeto para manipular memoria
   CtrlMemoria     cm;
   /// Puntero a resuelve sistema lineal
   ResuelveSistemaLineal *pRSL;

   /// Subdominio
   MatrizDensa subdominio;
   /// Particion del subdominio
   VectorInt   part;
   /// Toleracia usada en los calculos
   double      eps;
#ifdef DEPURAR
   /// Visualizaci�n de matrices y vectores generados
   bool visualiza;
#endif


   /// N�mero de nodos
   int num_nodos;
   /// N�mero de elementos
   int num_elementos;
   /// Numero de nodos desconocidos
   int num_nodos_interiores;
   /// N�mero de nodos en la frontera interior
   int num_nodos_frontera_interior;
   /// N�mero m�ximo de nodos en la frontera
   int num_maximo_nodos_en_frontera;
   /// Dimension del dominio
   int dimension;


   /// Vector soluci�n en los nodos interiores
   Vector          *pX;
   /// Vector soluci�n en todo el dominio
   VectorExtendido *pXto;
   /// Vector de carga local
   Vector          *pB;
   /// Vector de carga global
   VectorExtendido *pBA;
   /// Matriz de interiores-interiores o matriz de carga
   MatrizBandDisp  *pAII;
   /// Matriz de interiores-frontera interior
   MatrizBandDisp  *pAGI;
   /// Matriz de frontera interior-interiores
   MatrizBandDisp  *pAIG;
   /// Matriz de frontera interior-frontera interior
   MatrizBandDisp  *pAGG;



   /// Variables temporales usadas en Multiplica_Matriz_Schur
   Vector *py;
   Vector *py1;
   Vector *py2;
   Vector *pxb;
   Vector *psx;
   Vector *pBi;
   Vector *pSU;
   Vector *pUS;

   /// Nodos de frontera interior
   int    *nd_s;


   /// Matriz temporal para formar la matriz de rigidez
   MatrizDensa   *pXA;

   /// N�mero de puntos usados en la cuadratura Gaussiana
   int num_puntos_cuadratura;


   /// Genera el vector de carga
   void generaVectorCarga(void);

   /// Genera la matriz de interiores-interiores o matriz de carga
   void generaMatriz_AII(void);

   /// Genera la matriz de interiores-frontera interior
   void generaMatriz_AIG(void);

   /// Genera la matriz de frontera interior-interiores
   void generaMatriz_AGI(void);

   /// Genrea la matriz de frontera interior-frontera interior
   void generaMatriz_AGG(void);


   /// Inicializa el complemento de Schur
   void inicializaComplementoSchur(void);


public:

   /// Constructor de la clase
   DDM_SchurSubdominio(void)
   {
      pBA = NULL;
      pXto = NULL;
      pX = NULL;
      pB = NULL;
      pAII = NULL;
      pAGI = NULL;
      pAIG = NULL;
      pAGG = NULL;

      py = NULL;
      py1 = NULL;
      py2 = NULL;
      pxb = NULL;
      psx = NULL;
      pBi = NULL;
      pSU = NULL;
      pUS = NULL;

      pRSL = NULL;

      nd_s = NULL;

      num_puntos_cuadratura = 2;
      pXA = NULL;
#ifdef DEPURAR
      visualiza = false;
#endif
   }

   ~DDM_SchurSubdominio()
   {
      //~ liberaMemoria();
   }

   /// Destructor de la clase
   void liberaMemoria(void)
   {
      if (pBA) delete pBA;
      if (pX) delete pX;
      if (pB) delete pB;
      if (pXto) delete pXto;
      if (pAII)
      {
         pAII->liberaMemoria();
         delete pAII;
      }
      if (pAGI)
      {
         pAGI->liberaMemoria();
         delete pAGI;
      }
      if (pAIG)
      {
         pAIG->liberaMemoria();
         delete pAIG;
      }
      if (pAGG)
      {
         pAGG->liberaMemoria();
         delete pAGG;
      }

      if (py) delete py;
      if (py1) delete py1;
      if (py2) delete py2;
      if (pxb) delete pxb;
      if (psx) delete psx;
      if (pBi) delete pBi;
      if (pSU) delete pSU;
      if (pUS) delete pUS;

      if (pRSL) delete pRSL;

      if (nd_s) cm.liberaMemoriaINT(nd_s);

      if (pGe2D) delete pGe2D;
      pGe2D = NULL;

   }

   /// Retorna el numero de nodos desconocidos
   int retornaNumeroNodosDesconocidos(void)
   {
      return num_nodos_interiores;
   }

   /// Proyeccion de la frontera interior
   void proyeccionFronteraInterior(void);

   /// Soluci�n de los nodos interiores
   void solucionNodosInteriores(void);

   /// Retorna la solucion del subdominio
   void retornaSolucion(ldouble **coord, ldouble *val);

   /// Retorna error en el subdominio
   ldouble retornaErrorSubdominio(void);

   /// Asigna vector frontera interior
   void asignaVectorFronteraInterior(ldouble *val)
   {
      for (int i = 0; i < num_nodos_frontera_interior; i++) pUS->asigna(i,val[i]);
   }

   /// Retorna vector frontera interior
   void retornaVectorFronteraInterior(ldouble *val)
   {
      for (int i = 0; i < num_nodos_frontera_interior; i++) val[i] = pSU->retorna(i);
   }


   /// Retorna el contenido de la posicion i del vector Bi de Schur
   void contenidoVectorSchur(ldouble *val)
   {
      for (int i = 0; i < num_nodos_frontera_interior; i++) val[i] = pBi->retorna(i);
   }

   /// Retorna el contenido del vector de carga de la frontera interior
   void contenidoVectorCarga(ldouble *val)
   {
      int i, j;
      for (j = 0; j < num_nodos_frontera_interior; j++)
      {
         for (i = 0; i < num_nodos; i++) if (nd_s[i] == j ) val[j] = pBA->retorna(i);
      }
   }

   /// Puntero al problema a resolver
   void punteroProblema(Problema_2D *pr)
   {
      pPr2D = pr;
   }

   /// Retorna los nodos de la frontera del subdominio
   void retornaNodosFronteraSubdominio(ldouble **coord, int tam);

   /// Inicializa el subdominio
   void inicializaSubdominio(ldouble *sub, int *part, ldouble eps);

   /// Genera las matrices del m�todo
   void generaMatrices(void)
   {
      // Genera el vector de carga
      generaVectorCarga();
      // Genera la matriz de carga
      generaMatriz_AII();
      // Genera la metriz de frontera interior-interiores
      generaMatriz_AGI();
      // Genera la matriz de interiores-frontera interior
      generaMatriz_AIG();
      // Genera la matriz de frontera interior-frontera interior
      generaMatriz_AGG();

      // Borra matriz temporal
      delete pXA;
      pXA = NULL;

      // Libera en la geometria el soporte de nodos por elemento
      pGe2D->liberaSoporteNodos();
      // Inicializa el complemento de Schur
      inicializaComplementoSchur();
   }

   //~ /// Retorna el N�mero m�ximo de nodos de frontera
   //~ int retornaNumeroMaximoNodosFrontera(void)
   //~ { return (part.retorna(0) * 2 + part.retorna(1) * 2); }

#ifdef DEPURAR
   /// Activa/desactiva la visualizaci�n de matrices y vectores generados
   /** @param cad Puntero a la cadena de error */
   void activaVisualiza(bool st)
   {
      visualiza = st;
   }
#endif

};


#endif

