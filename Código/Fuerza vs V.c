/*---              Modelo de Tomlinson para Friccion                   ---*/
/*---               Una part�cula en una Dimensi�n                     ---*/
/*---                           16/04/2012                             ---*/
/*---        Trabaja con archivos de entrada con 20000 puntos          ---*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<conio.h>


FILE *textfi;
/* Colocar el nombre del archivo */
#define filename "c:/smanzi/2021/Tomlinson/Datos/FvsV_T298_cos1_E0317_k0.01_a029_f1E7(a).dat"
//#define filename "/home/smanzi/Tomlinson/promedios/E7000_T50.dat"

FILE *pun;
/*archivo lee*/
#define archi  "c:/smanzi/2021/Tomlinson/Matrices/Matrix_DE_cos1_E0317_k0.01_a029.dat"
const int DATOS=30003;

const double T=298;     /* Temperatura */
const double kB=1.3806503e-23;     /* Constante de Boltzmann */
const double conversion=1.60217646e-19;

const double k_real=0.01;   /* el valor de la constante del resorte en N/m */
//const double K=k_real*(1e-18)/(kB*T);  /* constante del cantilever*10^(-18)/(k_B T) */

//const double E0=0.140260976722*conversion/(kB*T);
const double a=0.29;    /* espaciamiento interatomico en nm */
const int saltos=5000; /* cantidad de saltos del tip */

/* cantidad de saltos del tip */
const double factor=1e7;   /* preexponencial en la probabilidad de salto */
const int cant_puntos=20;   //13
//const double pi=3.1415926535897932384626433832795;


int i,j,ii,jj,n,nprom,particula,lado,pos,intervalo,param,tip_0;
double random1,random2,DCM,xmedio,ZZ,DE,prob_der,prob_izq,prob_quieto,valor,ZZ1,R;
double tipinicial,tipfinal,valorMM,valorMM2,Kprog,tiempo,Dtiempo,pasoO,unidades;
double constante1,constante2,Ka,Delta,lugar,prob_total;
double pozo_i,tipfinal_i,tipfinal_d,pozo_d,p_i,auxiliar,coef_izq,coef_der,E_izq,E_der;
double distancia[saltos+1],distance,distancia0,DE_i,DE_d,Temp_rel;
int veces[saltos+1],saltosizquierda,chequeo,saltostotales;

double x1,x2,x1i,salto,velocidad;
double f1izq,f1der,f2izq,f2der,x1medio,x2medio,x1izq,x2izq,x1der,x2der,f1int,f2int,Fuerza_media;
int iiterar,xyz,fza_max,fm;

int ijk,tip;
double x1mediob,pend1,pend2,maximo,E02,K2,cte1;
double log_velocity, vel[cant_puntos];

double valores[cant_puntos],valores2[cant_puntos];
double cantilever[DATOS],mi[DATOS],MI[DATOS],centro[DATOS],MD[DATOS],md[DATOS];

double Dx,Dz,Z_inicial,Z_ref;
int saltoizquierda[DATOS];

int central,relativo,RR;
double porcentaje,tope_ZZ,centralita,DE1,E_referencia,limite_derecho,limite_izquierdo;

int nn;


void grabar(void);
void grabar1(void);
void calculos(void);
void ENTRADA(void);
double funcion(double punto);
double MinimoIzquierdo(void);
double MaximoIzquierdo(void);
double MinimoCentral(void);
double MaximoDerecho(void);
double MinimoDerecho(void);

/*---------------------------------------------------------------------*/
/*----        Generador  de  Nut  de  n�meros  aleatorios          ----*/

  const double long_max=2147483648.0;
  const double dos_long_max=4294967296.0;
  int  mas1_55[55], mas31_55[55], j_ran, i00, j00;
  long int  randi[55];
  void auxrnd(void);
  void randomizar(void);
  double ran01(void);
  unsigned int condi_55(int);
/*---------------------------------------------------------------------*/

void auxrnd(void)
  {   int rre;
     for(rre=0 ; rre<55 ; rre++)
     {mas1_55[rre]=condi_55(rre+1); mas31_55[rre]=condi_55(rre+31); }    }
void randomizar(void)
 {j_ran=0;
  for(i00=0;i00<55;i00++) randi[i00]=rand();
  for(i00=0;i00<1000;i00++)
  for(j00=0;j00<1000;j00++)
   {j_ran=mas1_55[j_ran]; randi[j_ran]=randi[j_ran]+randi[mas31_55[j_ran]];} }
double ran01(void)
   {  j_ran=mas1_55[j_ran];
      randi[j_ran]=randi[j_ran]+randi[mas31_55[j_ran]];
      return (((double)randi[j_ran]+long_max)/dos_long_max);   }
unsigned int condi_55(int z1)
  { int w;
if((z1>=0)*(z1<55)) w=z1; if (z1<0) w=55+z1; if (z1>=55)  w=z1-55; return w; }
/*---------------------------------------------------------------------*/


/*---------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void minimos_finos(void){

Z_ref=ZZ-(double)(tip)*Dz;

if(Z_ref>tope_ZZ){
      Z_ref=tope_ZZ;
      ZZ=Z_ref+(double)(tip)*Dz;
      prob_izq=0.0; prob_der=2.0;
      tipinicial=-centro[0]+(double)(tip)*Dx;
      Fuerza_media+=(ZZ-tipinicial); fm++;
      } /* sobrepaso de la fuerza limite  */
else{ /* 1 */

/*  Busca minimo de la particula  */
relativo=(int)((Z_ref-cantilever[0])/pasoO);
porcentaje=(Z_ref-cantilever[relativo])/pasoO;


tipinicial=MinimoCentral()+(double)(tip)*Dx;
Fuerza_media+=(ZZ-tipinicial);  fm++;

pozo_i=MinimoIzquierdo()+(double)(tip)*Dx;
DE_i=MaximoIzquierdo();
DE_d=MaximoDerecho();
pozo_d=MinimoDerecho()+(double)(tip)*Dx;


            if(Z_ref>limite_izquierdo){ prob_izq=0.0; }
            else{
                 prob_izq=exp(-DE_i);
                 }

            if(Z_ref<limite_derecho){ prob_der=0.0; }
            else{
                 prob_der=exp(-DE_d);
                 }
     }  /* cierra else 1 */

                         }
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void mayores(void){

minimos_finos();

prob_total=(prob_izq+prob_der);
p_i=factor*prob_total/R;
random1=ran01();

if(random1<p_i){ /* Intento saltar */
p_i=prob_izq/prob_total;
random2=ran01();
if(random2<p_i){ /* Salto a la izquierda aceptado */
                     distancia[tip]+=(ZZ-tipinicial);
                     veces[tip]++;
                     tip--;
                     saltosizquierda++;
                     saltostotales++;
                }
else{  /* Salto a la derecha aceptado */
                     distancia[tip]+=(ZZ-tipinicial);
                     veces[tip]++;
                     tip++;
                     saltostotales++;
      } /* Cierra else de intento */

                }  /* Cierra if de intento */

random1=1.0-ran01();
Dtiempo=-log(random1)/R;
ZZ+=velocidad*Dtiempo;

                       }
/*-------------------------------------------------------------------*/
double MinimoIzquierdo(void){

lugar=porcentaje*(mi[relativo+1]-mi[relativo])+mi[relativo];
return(lugar);
                                          }
/*-------------------------------------------------------------------*/
double MaximoIzquierdo(void){

lugar=porcentaje*(MI[relativo+1]-MI[relativo])+MI[relativo];
return(lugar);
                                          }
/*-------------------------------------------------------------------*/
double MinimoCentral(void){

lugar=porcentaje*(centro[relativo+1]-centro[relativo])+centro[relativo];
return(lugar);
                                        }
/*-------------------------------------------------------------------*/
double MaximoDerecho(void){

lugar=porcentaje*(MD[relativo+1]-MD[relativo])+MD[relativo];
return(lugar);
                                        }
/*-------------------------------------------------------------------*/
double MinimoDerecho(void){

lugar=porcentaje*(md[relativo+1]-md[relativo])+md[relativo];
return(lugar);
                                        }
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

void main(void){
/*randomize();*/
xyz=0;
ENTRADA();
//grabar1();

for(i=0;i<cant_puntos;i++){ log_velocity=2.0-0.5*i; vel[i]=1000.0*expl(log_velocity); }

printf("Temperatura = %f \n",T);

tope_ZZ=-cantilever[0];

pasoO=cantilever[3]-cantilever[2];

auxrnd(); randomizar();

i=0;
RR=1;
while(RR){
  if(md[i]){ limite_derecho=cantilever[i]; RR=0; }
  else{ i++; }
           }

i=0;
RR=1;
while(RR){
   if(mi[i]){ i++; }
   else{ limite_izquierdo=cantilever[i-1]; RR=0; }
          }



if(limite_izquierdo<0){ Z_inicial=0.0; }
else{ Z_inicial=limite_izquierdo; }

Dx=a; //MD[i-1]-centro[0];
Dz=a; //cantilever[i-2]-cantilever[0];
Dx=Dz;
/*------------------------------------*/
/*------------------------------------*/

R=2.0*factor;



tip_0=saltos-1;


ii=0;
for(xyz=0;xyz<cant_puntos;xyz++){

ENTRADA();
velocidad=vel[xyz];
printf("Velocidad = %f \n",log(vel[xyz]/1000.0));

Fuerza_media=0.0;  fm=0;
for(tip=0;tip<saltos;tip++){ distancia[tip]=0.0; veces[tip]=0.0; }

tiempo=0.0;
ZZ=Z_inicial;
Z_ref=ZZ;
relativo=(int)((Z_ref-centro[0])/pasoO);
porcentaje=(Z_ref-cantilever[relativo])/pasoO;


tip=0;
saltosizquierda=0;
saltostotales=0;

tipinicial=MinimoCentral()+(double)(tip)*Dx;
pozo_i=MinimoIzquierdo()+(double)(tip)*Dx;
DE_i=MaximoIzquierdo();
DE_d=MaximoDerecho();
pozo_d=MinimoDerecho()+(double)(tip)*Dx;


while(tip<saltos){ mayores(); }

distance=0.0;
valores2[ii]=0;
for(tip=1;tip<tip_0;tip++){
      distance+=(distancia[tip]/((double)(veces[tip]))) ;
      valores2[ii]+=veces[tip];
                                    }
valores[ii]=distance/((double)(saltos)-2.0);
valores2[ii]/=((double)(saltos)-2.0);
ii++;

grabar();

}

//grabar();
                }  /*- Cierra el main -*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void ENTRADA(void){



if((pun=fopen(archi,"r"))==NULL){ printf("error red"); getch();}
else{ /*1*/

for(i=0;i<DATOS;i++){
  	fscanf(pun,"%le %le %le %le %le %le \n",&cantilever[i],&mi[i],&MI[i],&centro[i],&MD[i],&md[i]);
         	     }
	fclose(pun);
   }/*1*/

  unidades=conversion/(kB*T);

for(i=0;i<DATOS;i++){
  	MI[i]=unidades*MI[i];
   MD[i]=unidades*MD[i];
         	     }


                   }
/*-------------------------------------------------------------------*/
void grabar1(void){
   if ( (textfi =fopen (filename,"wt"))== NULL ){
	      printf("Error al abrir el archivo\n");    }

               fprintf(textfi,"Velocidad Fza_Max Fza_Media Saltos_Izq \n");
       fclose(textfi);
                  }

/*-------------------------------------------------------------------*/
void grabar(void){
   if ( (textfi =fopen (filename,"a"))== NULL ){
	      printf("Error al abrir el archivo\n");    }
//   for(nn=0;nn<cant_puntos;nn++){
//               fprintf(textfi,"%f %f %f %f %f %f \n",ZZ,pozo_i,tipfinal_i,tipinicial,tipfinal_d,pozo_d);
               fprintf(textfi,"%f %1.12f %1.12f %1.12f \n",log(vel[xyz]/1000.0),k_real*valores[xyz],k_real*Fuerza_media/(double)(fm),(double)(saltosizquierda)/(double)(saltostotales)); //}
       fclose(textfi);
                  }


