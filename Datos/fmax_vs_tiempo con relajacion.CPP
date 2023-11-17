/*---              Modelo de Tomlinson para Friccion                   ---*/
/*---               Una partícula en una Dimensión                     ---*/
/*---                           16/04/2012                             ---*/
/*---        Trabaja con archivos de entrada con 20000 puntos          ---*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<conio.h>



FILE *pun;
/*archivo lee*/
#define archi  "c:/smanzi/2023/Multiple/Matrices/Matrix_DE_cos1_E0317_k1_a029.dat"


const int DATOS=20005;


const int N=1000;

const double Temperatura=1000.0;     /*  Temperatura */
const double kB=1.3806503e-23;     /* Constante de Boltzmann */
const double conversion=1.60217646e-19;

const double k_real=1.0;   /* el valor de la constante del resorte en N/m */
//const double K=k_real*(1e-18)/(kB*T);  /* constante del cantilever*10^(-18)/(k_B T) */

//const double E0=0.140260976722*conversion/(kB*T);
const double a=0.29;    /* espaciamiento interatomico en nm */
const int saltos=1; //10000; /* cantidad de saltos del tip */
const double Velocity=0.0;

const int ancho=0;     /* variación entre todos juntos =0 y todos bien separados uniformemente =1; para ancho=1 la máxima separación es a/2 */
const int relajando=5000; /* cantidad de relajaciones para los tip */

const int factor_grabar=1*N; /* cantidad de intentos entre cada grabación */


/* cantidad de saltos del tip */
const double factor=1e7;        /* preexponencial en la probabilidad de salto */
const int cant_puntos=1;   //13
//const double pi=3.1415926535897932384626433832795;

FILE *fptm;
char gfmax[50];


FILE *fpt;
char gfuerza[50];

FILE *fptd;
char gdistribucion[50];

int i,j,ii,jj,n,nprom,particula,lado,pos,intervalo,tip_0,saltitos,anchito,ir,irr,controlando,grabo;
double NN,ANCHO,random1,random2,DCM,xmedio,DE,prob_der[N],prob_izq[N],prob_der_s[N],prob_izq_s[N],prob_quieto,valor,ZZ1,R;
double tipinicial[N],tipfinal,valorMM,valorMM2,Kprog,tiempo,tiempo2,Dtiempo,pasoO,unidades,Dtiempo_salto;
double constante1,constante2,Ka,Delta,lugar,prob_total;
double pozo_i[N],tipfinal_i,tipfinal_d,pozo_d[N],p_i,auxiliar,coef_izq,coef_der,E_izq,E_der,Primer_Z[N];
double distance,distancia0,DE_i[N],DE_d[N],Temp_rel,velocidad_real,velocidad_real2,diferencia,dist_total;
int saltosizquierda[N],saltosderecha[N],chequeo,obligado,pasado;

/* saltos deberia ser saltos/N por 2 como factor de resguardo */

double x1,x2,x1i,salto,velocidad,SS;
double f1izq,f1der,f2izq,f2der,x1medio,x2medio,x1izq,x2izq,x1der,x2der,f1int,f2int;
int iiterar,xyz,fza_max,N1,ijk,gg;

int tip[N],Max_Desp[N];
double x1mediob,pend1,pend2,maximo,E02,K2,cte1,fuerza_total;
double log_velocity, vel[cant_puntos],ZZ[N];

double valores[cant_puntos],valores2[cant_puntos];
double cantilever[DATOS],mi[DATOS],MI[DATOS],centro[DATOS],MD[DATOS],md[DATOS];

double Dx,Dz,Z_inicial,Z_ref[N],dist_salto,T;


int central,relativo,RR;
double porcentaje,tope_ZZ,centralita,DE1,E_referencia,limite_derecho,limite_izquierdo;

int nn,exito,contacto,dd,dm;
double DDt;


void grabar_fmaxima(void);
void calculos(void);
void ENTRADA(void);
double funcion(double punto);
double MinimoIzquierdo(void);
double MaximoIzquierdo(void);
double MinimoCentral(void);
double MaximoDerecho(void);
double MinimoDerecho(void);
void mayores_Relajacion(void);
void minimos_finos_Relajacion(void);
void Relajacion(void);
void mide_fuerza(void);
void mide_distribucion(void);

/*---------------------------------------------------------------------*/
/*----        Generador  de  Nut  de  números  aleatorios          ----*/

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

      pasado=0;
      diferencia=0.0;

for(contacto=0;contacto<N;contacto++){

    Z_ref[contacto]=ZZ[contacto]-(double)(tip[contacto])*Dz;
    if((Z_ref[contacto]-tope_ZZ)>diferencia){  diferencia=(Z_ref[contacto]-tope_ZZ);
                                               obligado=contacto;
                                               pasado=1; }
                                      }

if(pasado){
    for(contacto=0;contacto<N;contacto++){
    ZZ[contacto]-=diferencia;
    Z_ref[contacto]-=diferencia;
                                      }

     tiempo2+=(Dtiempo-diferencia/velocidad);
     tipinicial[obligado]=-centro[0]+(double)(tip[obligado])*Dx;


           } /* sobrepaso de la fuerza limite  */
else{ /* 1 */

tiempo2+=Dtiempo;
for(contacto=0;contacto<N;contacto++){
/*  Busca minimo de la particula  */
relativo=(int)((Z_ref[contacto]-cantilever[0])/pasoO);
porcentaje=(Z_ref[contacto]-cantilever[relativo])/pasoO;


tipinicial[contacto]=MinimoCentral()+(double)(tip[contacto])*Dx;

//pozo_i[contacto]=MinimoIzquierdo()+(double)(tip[contacto])*Dx;
DE_i[contacto]=MaximoIzquierdo();
DE_d[contacto]=MaximoDerecho();
//pozo_d[contacto]=MinimoDerecho()+(double)(tip[contacto])*Dx;


            if(Z_ref[contacto]>limite_izquierdo){ prob_izq[contacto]=0.0; }
            else{
                 if(tip[contacto]){ prob_izq[contacto]=exp(-DE_i[contacto]); }
                 else{ prob_izq[contacto]=0.0; }
                 }

            if(Z_ref[contacto]<limite_derecho){ prob_der[contacto]=0.0; }
            else{
                 prob_der[contacto]=exp(-DE_d[contacto]);
                 }
     }  /* cierra else 1 */
                            } /* cierra for de contacto */
                         }
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void mayores(void){

minimos_finos();

if(pasado){
                     saltosderecha[obligado]++;
                     exito=0;
                     //printf("Tip = %d , \n",obligado);


                 for(ijk=0;ijk<N;ijk++){
                 Z_ref[ijk]=ZZ[ijk]-(double)(tip[ijk])*Dz;
                 relativo=(int)((Z_ref[ijk]-cantilever[0])/pasoO);
                 porcentaje=(Z_ref[ijk]-cantilever[relativo])/pasoO;
                 tipinicial[ijk]=MinimoCentral()+(double)(tip[ijk])*Dx;
                                       }

                     contacto=obligado;
                     grabar_fmaxima();
                     Dtiempo_salto=0.0;
                     tip[obligado]++;
                     saltitos++;

           }  /* cierra if de pasado = 1 */

else{

prob_total=(prob_izq[0]+prob_der[0]);

for(contacto=1;contacto<N;contacto++){

    prob_total+=(prob_izq[contacto]+prob_der[contacto]); }


p_i=factor*prob_total/R;
random1=ran01();

if(random1<p_i){ /* Intento saltar */

prob_der_s[0]=prob_der[0]/prob_total;
prob_izq_s[0]=prob_der_s[0]+prob_izq[0]/prob_total;


for(contacto=1;contacto<N;contacto++){
    prob_der_s[contacto]=prob_izq_s[contacto-1]+prob_der[contacto]/prob_total;
    prob_izq_s[contacto]=prob_der_s[contacto]+prob_izq[contacto]/prob_total;
                                      }


exito=1; contacto=0;
random1=ran01();
while(exito){


if(random1<prob_der_s[contacto]){  /* Salto a la derecha aceptado */

                     saltosderecha[contacto]++;
                     exito=0;
                     //printf("Tip que salta = %d \n",contacto);
                     if(tip[contacto]>Max_Desp[contacto]){
                        Max_Desp[contacto]=tip[contacto];
                        Dtiempo_salto=0.0;
                          }
                         grabar_fmaxima();
                         tip[contacto]++;
                     saltitos++;

                }
else{  if(random1<prob_izq_s[contacto]){/* Salto a la izquierda aceptado */
                     if(tip[contacto]){
                        tip[contacto]--;
                        saltosizquierda[contacto]++;
                        saltitos--;
                        grabar_fmaxima();
                                       }
                     exito=0;

            }
        else{
             contacto++;
             if(contacto==N){ exito=0; }}
      } /* Cierra else de intento */

                } /* cierra while de exito */


                }  /* Cierra if de intento */
       } /* cierra else de pasado=0 */



random1=1.0-ran01();
Dtiempo=-log(random1)/R;
tiempo+=Dtiempo;
Dtiempo_salto+=Dtiempo;
dist_salto=velocidad*Dtiempo;
for(contacto=0;contacto<N;contacto++){
    ZZ[contacto]+=dist_salto;
                                      }

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
randomize();
N1=N-1;
xyz=0;

NN=(double)(N);
ANCHO=(double)(ancho);

T=Temperatura;
ENTRADA();
tiempo=0.0;
tiempo2=0.0;
Dtiempo_salto=0.0;
gg=0;
anchito=ancho;

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

R=2.0*factor*(double)(N);

tip_0=saltos-1;

ii=0;
for(xyz=0;xyz<cant_puntos;xyz++){


velocidad=1000.0*expl(Velocity);
printf("Velocidad = %f \n",velocidad);


tiempo=0.0;


dd=1; dm=0;

/*-------------------------------------------------------------------------*/
/* Aca se eligen las posiciones iniciales de los cantilever aleatoriamente */


if(anchito){ /* Para aleatorios */


SS=0.5*a;
if(SS>(-cantilever[0])){ SS=-cantilever[0]; } /* sin esta restricción, cuando el resorte es muy duro no encontrará el mínimo en algunas ocasiones */
//SS=-cantilever[0];
for(ijk=0;ijk<N;ijk++){
    random1=ran01();
    Z_ref[ijk]=SS*(2.0*random1-1.0);
                        }
             }
else{  /* Todos juntos */
    for(ijk=0;ijk<N;ijk++){   Z_ref[ijk]=0.0; }
     }
   irr=-1;

/*-------------------------------------------------------------------------*/

for(contacto=0;contacto<N;contacto++){
    tip[contacto]=0;
    Max_Desp[contacto]=0;
                                      }


for(contacto=0;contacto<N;contacto++){

relativo=(int)((Z_ref[contacto]-cantilever[0])/pasoO);
porcentaje=(Z_ref[contacto]-cantilever[relativo])/pasoO;

tipinicial[contacto]=MinimoCentral();
pozo_i[contacto]=MinimoIzquierdo();
DE_i[contacto]=MaximoIzquierdo();
DE_d[contacto]=MaximoDerecho();
pozo_d[contacto]=MinimoDerecho();
                                          }
saltitos=0;

controlando=1;
printf("Comienza Relajacion \n");
//mide_distribucion();
Relajacion();
//mide_distribucion();

/*-------------------------------------------------------------------------*/

T=Temperatura;
ENTRADA();

printf("Comienza Friccion \n");

for(contacto=0;contacto<N;contacto++){
    ZZ[contacto]=Z_ref[contacto];
    tip[contacto]=0;
    Max_Desp[contacto]=0;
    saltosizquierda[contacto]=0;
    saltosderecha[contacto]=0; }


for(contacto=0;contacto<N;contacto++){

relativo=(int)((ZZ[contacto]-cantilever[0])/pasoO);
porcentaje=(ZZ[contacto]-cantilever[relativo])/pasoO;

tipinicial[contacto]=MinimoCentral();
pozo_i[contacto]=MinimoIzquierdo();
DE_i[contacto]=MaximoIzquierdo();
DE_d[contacto]=MaximoDerecho();
pozo_d[contacto]=MinimoDerecho();
                                          }
saltitos=0;


  grabar_fmaxima();
  grabo=0;

while(saltitos<saltos){
  //     if(grabo>factor_grabar){   grabar_fmaxima(); grabo=0; }
                 grabo++;
                 mayores(); }

//grabar();

}

                }  /*- Cierra el main -*/

/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
void Relajacion(void){


for(irr=0;irr<relajando;irr++){
for(ir=0;ir<N;ir++){
  mayores_Relajacion();
                     }
                     }
   mide_fuerza();
   mide_distribucion();
 //if(controlando==irr){  mide_distribucion(); controlando*=10; }

                      }
/*-------------------------------------------------------------------*/
void mayores_Relajacion(void){


ijk=(int)((double)(N)*ran01());

if(Z_ref[ijk]>tope_ZZ){ prob_der[ijk]=1.0; prob_izq[ijk]=0.0; Z_ref[ijk]=tope_ZZ; }
else{
    if(Z_ref[ijk]<(-tope_ZZ)){ prob_der[ijk]=0.0; prob_izq[ijk]=1.0; Z_ref[ijk]=-tope_ZZ; }
    else{
         relativo=(int)((Z_ref[ijk]-cantilever[0])/pasoO);
         porcentaje=(Z_ref[ijk]-cantilever[relativo])/pasoO;
         tipinicial[ijk]=MinimoCentral();
         DE_i[ijk]=MaximoIzquierdo();
         DE_d[ijk]=MaximoDerecho();
         if(DE_i[ijk]){ prob_izq[ijk]=exp(-DE_i[ijk]); }
         else{ prob_izq[ijk]=0.0; }
         if(DE_d[ijk]){ prob_der[ijk]=exp(-DE_d[ijk]); }
         else{ prob_der[ijk]=0.0; }
         }
     }

prob_total=(prob_izq[ijk]+prob_der[ijk]);
if(prob_total>1.0){ prob_izq[ijk]/=prob_total;
                    prob_der[ijk]/=prob_total; }

random1=ran01();
if(random1<prob_izq[ijk]){
                     Z_ref[ijk]+=a;
                          }
else{
     if(random1<(prob_izq[ijk]+prob_der[ijk])){
                     Z_ref[ijk]-=a;
                                               }
     }


                              }
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
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

/*-------------------------------------------------------------------*/
void grabar_fmaxima(void){

fuerza_total=0.0;
dist_total=0.0;


if(gg){
fuerza_total=0.0;
for(ijk=0;ijk<N;ijk++){
                 Z_ref[ijk]=ZZ[ijk]-(double)(tip[ijk])*Dz;
                 relativo=(int)((Z_ref[ijk]-cantilever[0])/pasoO);
                 porcentaje=(Z_ref[ijk]-cantilever[relativo])/pasoO;
                 tipinicial[ijk]=MinimoCentral()+(double)(tip[ijk])*Dz;

                 fuerza_total+=(ZZ[ijk]-tipinicial[ijk]); }

                 fuerza_total=k_real*fuerza_total/(double)(N);



 sprintf(gfmax,"c:/smanzi/2023/Multiple/Datos/Fmax_vs_t_Relajado_Ancho%1.0f_K%1.2f_Ntips%1.0f_Vel%1.2f_T%1.2f.dat",ANCHO,k_real,NN,Velocity,Temperatura);
 fptm=fopen(gfmax,"a");


       fprintf(fptm,"%1.12f %1.12f %d %1.12f \n",tiempo2,fuerza_total,saltitos,Dtiempo_salto);
       fclose(fptm);

       }

else{
   gg=1;
 sprintf(gfmax,"c:/smanzi/2023/Multiple/Datos/Fmax_vs_t_Relajado_Ancho%1.0f_K%1.2f_Ntips%1.0f_Vel%1.2f_T%1.2f.dat",ANCHO,k_real,NN,Velocity,Temperatura);
 fptm=fopen(gfmax,"wt");
               fprintf(fptm,"TiempoReal FuerzaFriccion Nro_salto Dt_salto \n");
       fclose(fptm);

     }

                  }

/*-------------------------------------------------------------------*/
void mide_fuerza(void){


fuerza_total=0.0;
for(ijk=0;ijk<N;ijk++){
                 relativo=(int)((Z_ref[ijk]-cantilever[0])/pasoO);
                 porcentaje=(Z_ref[ijk]-cantilever[relativo])/pasoO;
                 tipinicial[ijk]=MinimoCentral()+(double)(tip[ijk])*Dz;

                 fuerza_total+=(Z_ref[ijk]-tipinicial[ijk]); }

                 fuerza_total=k_real*fuerza_total/(double)(N);


sprintf(gfuerza,"c:/smanzi/2023/Multiple/Datos/Fuerza_Relajando_Ancho%d_K%f_T%f.dat",ancho,k_real,Temperatura);
fpt=fopen(gfuerza,"a");

fprintf(fpt,"%d %1.12f \n",irr+1,fuerza_total);

fclose(fpt);
                       }
/*-------------------------------------------------------------------*/
void mide_distribucion(void){

sprintf(gdistribucion,"c:/smanzi/2023/Multiple/Datos/Dist_Relax_Ancho%d__Rel%d_K%1.2f_T%f.dat",ancho,k_real,Temperatura,relajando);
fptd=fopen(gdistribucion,"w");

for(ijk=0;ijk<N;ijk++){
                 relativo=(int)((Z_ref[ijk]-cantilever[0])/pasoO);
                 porcentaje=(Z_ref[ijk]-cantilever[relativo])/pasoO;
                 tipinicial[ijk]=MinimoCentral();

fprintf(fptd,"%d %1.12f %1.12f \n",ijk,Z_ref[ijk],tipinicial[ijk]);
//fprintf(fptd,"%1.12f \n",Z_ref[ijk]);
                        }
fclose(fptd);

                             }
/*-------------------------------------------------------------------*/


