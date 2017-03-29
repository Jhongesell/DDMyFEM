// Programa Maestro - Esclavo

// Compilar usando
//   mpiCC.mpich -O1 me.cpp -o me -lm

// Correr usando 8 procesadores por ejemplo
//   mpirun.mpich -np 8 me




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>

// Funcion de Arnold
double F(double T,double A,double B);
//  Calcula  la  resonancia  de  algun  punto del espacio  de parametros
int Resonancias(int MaxPer,int Ciclos,int Trans,double A,double B, double Eps,int *Pr,int *Qr);



//#define VIS


// Programa Maestro-Esclavo
int main(int argc, char *argv[]) 
{
   int ME_id,MP_np;
   int ME_P, ME_L, ME_sw;
   
   MPI::Init(argc,argv);
   ME_id = MPI::COMM_WORLD.Get_rank();
   MP_np = MPI::COMM_WORLD.Get_size();

	// Revisa que pueda arrancar el esquema M-E
	if (MP_np < 2) 
   {
      printf("Se necesitan almenos dos procesadores para el esquema M-E\n");
      return 1;
   }
   
   int i;
   // Controlador del esquema M-E
	if (ME_id == 0) 
   {	
      // Maestro
		printf("Maestro-Esclavo, Numero de Esclavos %d\n",MP_np-1);

      long double *a = new long double [10];
      for(i = 0; i < 10; i++) a[i] = 10.0-i;
         
      // Aviso de envio de una nueva tarea al nodo L
      for(i = 1; i < MP_np; i++) MPI::COMM_WORLD.Send(a, 10, MPI::LONG_DOUBLE, i,1);

      
	} else {
      // Esclavos
   
      long double *b = new long double[10];
      // Recibo aviso de envio de una nueva tarea 
      MPI::COMM_WORLD.Recv(b, 10, MPI::LONG_DOUBLE, 0,1);
      for(i = 0; i < 10; i++) printf("%Lf\n",b[i]);

	}
	
	MPI::Finalize();
	return 0;
}


