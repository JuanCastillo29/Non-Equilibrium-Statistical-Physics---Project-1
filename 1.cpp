#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

using namespace std;


unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran,ig1,ig2,ig3;
#define NormRANu (2.3283063671E-10F)
#define Pi 3.14159265

#define NRan 10000
#define Tfin 100
#define F 0.01
#define dt 0.1
#define Gamma 100
#define NPart 1000
#define L 1


const double Phi=sqrt(2*Gamma*dt);

//Generador de números de Parisi-Rapuano.
void ini_ran(int SEMILLA)
{
    int INI,FACTOR,SUM,i;
    srand(SEMILLA);
    INI=SEMILLA;
    FACTOR=67397;
    SUM=7364893;
    for(i=0;i<256;i++)
    {
        INI=(INI*FACTOR+SUM);
        irr[i]=INI;
    }
    ind_ran=ig1=ig2=ig3=0;
}

float Random(void)
{
    float r;
    ig1=ind_ran-24;
    ig2=ind_ran-55;
    ig3=ind_ran-61;
    irr[ind_ran]=irr[ig1]+irr[ig2];
    ir1=(irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;
    return r;
}

//Generador de los números aleatorios según distribución gaussiana por el método de Box-Muller

double DistrGauss(void){
    double d1, d2;
    //Se generan los numeros planos necesarios
    d1 = Random();
    d2 = Random();
    //Se usa la fórmula para crear dos números independientes según una gaussiana
    return -sqrt(-2*log(d1))*cos(2*Pi*d2);
}


//En esta funcion creo las condiciones periodicas.
//r guarda la posicion de la particula en la caja definida y cont cuenta el número de vueltas dado para poder desplegar las trayectorias luego
void Contorno(double (&r)[NPart][2], int (&cont)[NPart][2]){
    for(int i=0; i<NPart; i++){
        for(int j=0; j<2; j++){
            if(r[i][j]>=L){
                cont[i][j] = cont[i][j]+1;
                r[i][j] = r[i][j]-L;
            }
            if(r[i][j]<0){
                cont[i][j] = cont[i][j]-1;
                r[i][j] = r[i][j]+L;
            }
        }
    }
}

void ActPosition(double (&r)[NPart][2], int (&cont)[NPart][2], double (&force)[NPart]){
    for(int i=0; i<NPart; i++){
        r[i][0] = r[i][0] + force[i]*dt + Phi*DistrGauss();
        r[i][1]= r[i][1] +Phi*DistrGauss();
    }
    Contorno(r, cont);
}

void Guardardatos(double (&r)[NPart][2], int t, int (&cont)[NPart][2]){
    char filename[20]; // Declaras un array de caracteres para almacenar el nombre del fichero
    // Formateas el nombre del fichero con sprintf
    sprintf(filename, "Evolucion/t%d.txt", t);

    FILE *fv1 = fopen(filename, "wt");
    if(fv1 == NULL)
            printf("No se pudo abrir el fichero");
    for (int n=0; n<NPart; n++)
        fprintf(fv1, "%lf \t %lf \n", r[n][0]+ L*cont[n][0], r[n][1]+ L*cont[n][1]);
    fclose(fv1);
}

void GuardarTrayectoria(double (&r)[NPart][2], int (&Contorn)[NPart][2], int t){
char filename[20], modo[20];
for(int i=0; i<NPart; i++){
  sprintf(filename, "Trayectorias/Particula %d.txt", i);
  if(t==0)
    sprintf(modo, "wt");
  else{
    sprintf(modo, "at");
  }
  FILE *f=fopen(filename, modo);
  if(f!=NULL){
    fprintf(f, "%lf\t%lf\n", r[i][0]+L*Contorn[i][0],r[i][1]+L*Contorn[i][1] );
    fclose(f);
  }
}
}

void Distribucion(double (&r)[NPart][2], int (&Cont)[NPart][2], int n){
char filename[20], filename2[20], filename3[20];
sprintf(filename, "Evolucion/t%d.txt", n);
sprintf(filename2, "Distribucion/t%d.txt", n);
sprintf(filename3, "Response function/t%d.txt", n);
FILE *f = fopen(filename, "rt"), *f0 = fopen("Evolucion/t0.txt", "rt"), *fw = fopen(filename2, "wt"), *fw2 = fopen(filename3, "wt");
if(f!=NULL && f0!=NULL && fw!=NULL && fw2!=NULL){
    double r0[NPart][2], rn[NPart][2];
    for(int i=0; i< NPart; i++){
        fscanf(f0, "%lf\t%lf\n", &r0[i][0],&r0[i][1]);
        fscanf(f, "%lf\t%lf\n", &rn[i][0],&rn[i][1]);
        fprintf(fw,"%lf\t%lf\n", rn[i][0] - r0[i][0], rn[i][1] - r0[i][1]);
        fprintf(fw2,"%lf\n", rn[i][0] - r0[i][0]);
    }
    fclose(f);
    fclose(f0);
    fclose(fw);
    fclose(fw2);
}
else{
    printf("No se pudo abrir el fichero.");
    if(f==NULL)
        printf("El primero.\n");
    if(f0==NULL)
        printf("El segundo\n");
    if(fw==NULL)
        printf("El tercero\n");
    if(fw2==NULL)
        printf("El cuarto\n");
}

}

void Trayectoria(void){
    double r[NPart][2], force[NPart];
    int Contorno[NPart][2];
    for(int i=0; i< NPart; i++){
        r[i][0]=Random();
        r[i][1]=Random();
        Contorno[i][0] = 0;
        Contorno[i][0] = 0;
        if(Random()<0.5)
            force[i]=-F;
        else
            force[i]=F;
    }
    for(int n =0; n<=Tfin/dt; n++){
        Guardardatos(r, n, Contorno);
        Distribucion(r, Contorno, n);
       //GuardarTrayectoria(r, Contorno, n);
        ActPosition(r, Contorno, force);
    }
}



void DiffConst(void){
char filename[20], filename2[20];
sprintf(filename, "Distancias/Gamma = %lf.txt", Gamma);
sprintf(filename2, "Distribucion/t%d.txt", int(Tfin/dt));
FILE *fdist =fopen(filename, "wt"), *ffin=fopen(filename2, "rt");
if(fdist == NULL || ffin ==NULL)
    printf("No se pudo abrir alguno de los ficheros");
else{
    double dx, dy, dist;
    for(int i=0; i<NPart; i++){
        fscanf(ffin, "%lf \t %lf\n", &dx, &dy);
        dist = (dx)*(dx) + (dy)*(dy);
        fprintf(fdist, "%lf\n", dist/(4*Tfin));
    }
    fclose(fdist);
    fclose(ffin);
}
}

int main(){
ini_ran(time(NULL));
Trayectoria();
//DiffConst();
return 0;
}
