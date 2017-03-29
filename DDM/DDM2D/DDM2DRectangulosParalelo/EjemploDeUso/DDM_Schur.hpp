//////////////////////////////////////////////////////////////////////////////////////////////
// Clase para definir los m�todos del m�todo de descomposici�n de dominio de subestructura- //
// ci�n Schur                                                                               //
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


#ifndef __DDM_Schur__
#define __DDM_Schur__


#include "Problema_2D.hpp"
#include "Geometria.hpp"
#include "DDM_SchurSubdominio.hpp"
#include "CtrlMemoria.hpp"
#include "Geometria_2DRectangulos.hpp"
#include "BCGM.hpp"


/// Clase para definir el m�todo descomposici�n de dominio de subestructuraci�n (Complemento de Schur)
/** @author Antonio Carrillo Ledesma
    @date verano 2009
    @version 1.0.0
    @bug No hay errores conocidos
*/
class DDM_Schur: public BCGM, public  MultOp, public ProductoPunto
{

private:


   // Configuraci�n necesaria para el m�todo CGM

   /// Producto punto, para el m�todo de CGM
   inline double productoPunto(Vector *u, Vector *v)
   {
      return u->productoPunto(v);
   }

   /// Multiplica Au=v, para el m�todo de CGM
   inline void multiplica(Vector *u, Vector *v)
   {
      int i, n, k, j;
      v->inicializa(0.0);
      // Calculo del vector mediante la multiplicaci�n de x por matriz de Schur Si de cada subdominio
      for (i = 0; i < num_sub; i++)
      {
         for (n = 0, j = 0; j < num_nodos_frontera_interior; j++)
         {
            for (k = 0; k < num_nodos_elemento; k++)
            {
               if (subdominio[j][k] == i)
               {
                  pXU[num_nodo_local[j][k]] = u->retorna(j);
                  n++;
               }
            }
         }
         funPar05(i, pXU, num_nodos_frontera_interior);

         for (j = 0; j < num_nodos_frontera_interior; j++)
         {
            for (k = 0; k < num_nodos_elemento; k++)
            {
               if (subdominio[j][k] == i) v->asigna(j,v->retorna(j) + pXU[num_nodo_local[j][k]]);
            }
         }
      }
   }

   /// Tama�o del sistema lineal, para el m�todo de CGM
   inline int tamano(void)
   {
      return num_nodos_frontera_interior;
   }



protected:

   /// Puntero a un arreglo de subdominios
   DDM_SchurSubdominio *subdom;
   /// Puntero a la definicion de la geometria
   Geometria_2D       *pGe2D;
   /// Puntero a la definicion del problema
   Problema_2D        *pPr2D;
   /// Objeto para manipular memoria
   CtrlMemoria cm;

   /// Particion
   VectorInt part;
   /// Tolerancia
   ldouble eps;
   /// N�mero de subdominios en la partici�n
   int num_sub;
   /// N�mero de nodos interiores
   int num_nodos_interiores;
   /// N�mero de nodos de frontera interior
   int num_nodos_frontera_interior;
   /// Dimensi�n de la descomposicion
   int dimension;

   /// N�mero de nodos por elemento
   int num_nodos_elemento;
   /// Valores soluci�n de la frontera interior
   Vector  *pXX;
   /// Valores intermedios de la solucion de la frontera interior
   ldouble  *pXU;
   /// Coordenadas de la frontera interior
   ldouble **coord_FI;
   /// Indica el numero de sumdominios que comparten el nodo
   int *multiplicidad;
   /// Indica al subdominio al que pertenece el nodo
   int **subdominio;
   /// Indica el n�mero de nodo local de los nodos de frontera interior del subdominio
   int **num_nodo_local;


   // Comportamientos a paralelizar

   /// Inicializa los subdominios
   virtual void funPar00(int ns)
   {
      subdom = new DDM_SchurSubdominio[ns];
      if (!subdom) cm.errorSolicitudMemoria("No hay memoria para los subdominios");
   }

   /// Recibe los datos para inicializar el subdominio
   virtual int funPar01(int ind, Problema_2D *pr, ldouble *sub, int *par, ldouble ep)
   {
      // Puntero al problema
      subdom[ind].punteroProblema(pr);
      // Inicializa el subdominio
      subdom[ind].inicializaSubdominio(sub, par, ep);
      // Actualiza el numero de nodos interiores desconocidos
      return subdom[ind].retornaNumeroNodosDesconocidos();
   }

   // Envia las coordenadas de la frontera interior
   virtual void funPar02(int ind, ldouble **cd, int tm)
   {
      subdom[ind].retornaNodosFronteraSubdominio(cd,tm);
   }

   // Genera las matrices
   virtual void funPar03(int ind)
   {
      subdom[ind].generaMatrices();
   }

   // Envia el vector de carga
   virtual void funPar04(int ind, ldouble *val1, ldouble *val2, int tm)
   {
      subdom[ind].contenidoVectorCarga(val1);
      subdom[ind].contenidoVectorSchur(val2);
   }

   // Proyeccion
   virtual void funPar05(int ind, ldouble *val, int tm)
   {
      subdom[ind].asignaVectorFronteraInterior(val);
      subdom[ind].proyeccionFronteraInterior();
      subdom[ind].retornaVectorFronteraInterior(val);
   }

   // Resuelve nodos interiores
   virtual void funPar06(int ind, ldouble *val, int tm)
   {
      subdom[ind].asignaVectorFronteraInterior(val);
      subdom[ind].solucionNodosInteriores();
   }

   // Retorna la solucion del subdominio
   virtual void funPar07(int ind, ldouble **cd, ldouble *vl, int tm)
   {
      subdom[ind].retornaSolucion(cd, vl);
   }

   virtual ldouble funPar08(int ind)
   {
      return subdom[ind].retornaErrorSubdominio();
   }

   virtual void funPar09(int ind)
   {
      subdom[ind].liberaMemoria();
   }

   /// Inicializa los subdominios
   void inicializaSubdominios(void);

   /// Conoce los nodos de la frontera interior
   void conoceFronteraInterior(void);

   /// Resuelve subdominios
   void resuelveSubdominios(void);

   /// Resuelve nodos interiores
   void resuelveNodosInteriores(void);


public:

   DDM_Schur(void) : BCGM(*(MultOp*) this, *(ProductoPunto*) this, 1000, 1e-5)
   {
      subdom = NULL;
      pXX = NULL;
      pXU = NULL;
      coord_FI = NULL;
      multiplicidad = NULL;
      subdominio = NULL;
      num_nodo_local = NULL;
      num_sub = 0;
   }

   /// Destructor de la clase
   ~DDM_Schur(void)
   {
      if (subdom)
      {
         int i;
         for (i = 0; i < num_sub; i++) funPar09(i);
         delete []subdom;
         subdom = NULL;
      }
      if (pXX) delete pXX;
      if (pXU) cm.liberaMemoriaLDOUBLE(pXU);
      if (coord_FI) cm.liberaMemoriaLDOUBLE(num_nodos_frontera_interior,coord_FI);
      if (multiplicidad) cm.liberaMemoriaINT(multiplicidad);
      if (subdominio) cm.liberaMemoriaINT(num_nodos_frontera_interior,subdominio);
      if (num_nodo_local) cm.liberaMemoriaINT(num_nodos_frontera_interior,num_nodo_local);
   }

   /// Constructor de la clase
   void inicializa(Problema_2D *pr, Geometria_2DRectangulos *ge, int mx, int my, double eps=1e-5)
   {
      if (mx < 2 || my < 2)
      {
         printf("\n\nError no se esta haciendo una descomposci�n fina valida (%d x %d)\n",mx,my);
         exit(1);
      }

      pPr2D = pr;
      pGe2D = ge;

      dimension = pGe2D->retornaDimension();
      part.redimensiona(dimension);
      part.asigna(0,mx);
      part.asigna(1,my);

      this->eps = eps;
   }


   /// Resuelve mediante el m�todo de descomposici�n de dominio
   void resuelveDDM(void)
   {
      // Inicializa subdominios
      inicializaSubdominios();
      // Conoce la frontera interior
      conoceFronteraInterior();
      // Resuelve subdominios
      resuelveSubdominios();
      // Resuelve nodos interiore
      resuelveNodosInteriores();
   }

   /// Retorna la soluci�n de cada punto de la discretizaci�n del subdominio
   /** @param sub N�mero del subdominio
       @param coord Coordenada
       @param val Valor
       @return N�mero de coordenadas y valores regresados */
   int retornaSolucion(int sub, ldouble **coord, ldouble *val);

   /// Retorna el m�ximo error de la soluci�n vs la soluci�n an�litica en el subdominio
   /** @param sub N�mero del subdominio
       @return El m�ximo error encontrado en todo el dominio*/
   ldouble retornaErrorSubdominio(int sub);

   /// Visualiza la solicion del problema
   void visualizaSolucion(void);

   /// Graba la soluci�n a un archivo
   /** @param arch Nombre del archivo en el cual se grabara la solucion del problema */
   void grabaSolucion(const char *arch);

   /// Retorna el m�ximo error de la soluci�n vs la soluci�n an�litica
   /** @return El m�ximo error encontrado en todo el dominio */
   ldouble error(void);


};

#endif
