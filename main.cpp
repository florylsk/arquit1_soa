#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

//realizado para 4000 objetos 
//debido a que c++ no admite arrays de longitud variable
struct Particulas{
    double posicionesX[4000];
    double posicionesY[4000];
    double posicionesZ[4000];
    double velocidadesX[4000];
    double velocidadesY[4000];
    double velocidadesZ[4000];
    double fuerzasX[4000];
    double fuerzasY[4000];
    double fuerzasZ[4000];
    double masas[4000];

};

const double G= 6.674*(pow(10,-11));

//variables globales
int num_objects;
int num_iteration;
int random_seed;
double size_enclosure;
double time_step;
double norma;
double accX;
double accY;
double accZ;


//crea el archivo init_config.txt
int init_config(Particulas &particulas){
    fstream file;
    file.open("init_config.txt",fstream::in | fstream::out | fstream::trunc);
    file<<fixed<<showpoint;
    file<<setprecision(3);
    file << size_enclosure<<" "<<time_step<<" "<<num_objects<<endl;
    for(int i = 0; i < num_objects; i++){

        file<<particulas.posicionesX[i]<<" "<<particulas.posicionesY[i]<<" "<<particulas.posicionesZ[i]<<" ";
        file<<particulas.velocidadesX[i]<<" "<<particulas.velocidadesY[i]<<" "<<particulas.velocidadesZ[i]<<" ";
        file<<particulas.masas[i]<<endl;

    }
    file.close();
    return 0;
}
//si 2 particulas han colisionado
void colision_particulas(Particulas &A,int posA, int posB){
    A.masas[posA] = A.masas[posA] + A.masas[posB];
    A.velocidadesX[posA] = A.velocidadesX[posA] + A.velocidadesX[posB];
    A.velocidadesY[posA] = A.velocidadesY[posA] + A.velocidadesY[posB];
    A.velocidadesZ[posA] = A.velocidadesZ[posA] + A.velocidadesZ[posB];
    A.posicionesX[posB]=0;
    A.posicionesY[posB]=0;
    A.posicionesZ[posB]=0;
    A.velocidadesX[posB]=0;
    A.velocidadesY[posB]=0;
    A.velocidadesZ[posB]=0;
    A.fuerzasX[posB]=0;
    A.fuerzasY[posB]=0;
    A.fuerzasZ[posB]=0;
    A.masas[posB]=0;
}
//generar las particulas
int generar_particulas(Particulas &particulas,mt19937_64 &generator,uniform_real_distribution<double> &dis,normal_distribution<double> &d){
    for(int i = 0; i < num_objects ; i++){
        particulas.posicionesX[i] = dis(generator);
        particulas.posicionesY[i] = dis(generator);
        particulas.posicionesZ[i] = dis(generator);
        particulas.velocidadesX[i]=0;
        particulas.velocidadesY[i]=0;
        particulas.velocidadesZ[i]=0;
        particulas.fuerzasX[i]=1;
        particulas.fuerzasY[i]=1;
        particulas.fuerzasZ[i]=1;
        particulas.masas[i] = d(generator);
    }

    return 0;

}
//comprobar si 2 particulas han colisionado
int comprobar_colision(Particulas &particulas, int pos){
    for(int i = 0; i < num_objects; i++){
        if(particulas.masas[i] != particulas.masas[pos]){
            if((particulas.posicionesX[i]-particulas.posicionesX[pos]) < 1 && (particulas.posicionesX[i]-particulas.posicionesX[pos]) > -1){
                if((particulas.posicionesY[i]-particulas.posicionesY[pos]) < 1 && (particulas.posicionesY[i]-particulas.posicionesY[pos]) > -1){
                    if((particulas.posicionesZ[i]-particulas.posicionesZ[pos]) < 1 && (particulas.posicionesZ[i]-particulas.posicionesZ[pos]) > -1){
                        return i;
                    }
                }
            }
        }
    }
    return -1;
}
//cada particula recorre toda la lista de particulas
//y se le suma la fuerza que le aplica cada una
void fuerza_gravitatoria(Particulas &particulas, int posA){
    for(int i=0;i<num_objects;i++){
        if(particulas.masas[posA] != particulas.masas[i]){
            if((particulas.posicionesX[posA]-particulas.posicionesX[i]) == 0){
                continue;
            }
            else if((particulas.posicionesY[posA]-particulas.posicionesY[i]) == 0){
                continue;
            }
            else if ((particulas.posicionesZ[posA]-particulas.posicionesZ[i]) == 0){
                continue;
            }
            else {
                norma=((std::sqrt(((particulas.posicionesX[i] - particulas.posicionesX[posA])*(particulas.posicionesX[i] - particulas.posicionesX[posA]))+((particulas.posicionesY[i] - particulas.posicionesY[posA])*(particulas.posicionesY[i] - particulas.posicionesY[posA]))+((particulas.posicionesZ[i] - particulas.posicionesZ[posA])*(particulas.posicionesZ[i] - particulas.posicionesZ[posA]))))*(std::sqrt(((particulas.posicionesX[i] - particulas.posicionesX[posA])*(particulas.posicionesX[i] - particulas.posicionesX[posA]))+((particulas.posicionesY[i] - particulas.posicionesY[posA])*(particulas.posicionesY[i] - particulas.posicionesY[posA]))+((particulas.posicionesZ[i] - particulas.posicionesZ[posA])*(particulas.posicionesZ[i] - particulas.posicionesZ[posA]))))*(std::sqrt(((particulas.posicionesX[i] - particulas.posicionesX[posA])*(particulas.posicionesX[i] - particulas.posicionesX[posA]))+((particulas.posicionesY[i] - particulas.posicionesY[posA])*(particulas.posicionesY[i] - particulas.posicionesY[posA]))+((particulas.posicionesZ[i] - particulas.posicionesZ[posA])*(particulas.posicionesZ[i] - particulas.posicionesZ[posA])))));
                particulas.fuerzasX[posA] += (G * particulas.masas[posA] * particulas.masas[i] * (particulas.posicionesX[i] - particulas.posicionesX[posA])) /
                               norma;
                particulas.fuerzasY[posA] += (G * particulas.masas[posA] * particulas.masas[i] * (particulas.posicionesY[i] - particulas.posicionesY[posA])) /
                               norma;
                particulas.fuerzasZ[posA] += (G * particulas.masas[posA] * particulas.masas[i] * (particulas.posicionesZ[i] - particulas.posicionesZ[posA])) /
                               norma;
            }
        }
    }
}

//a=f/m,v=v0+a*dt
void aceleracion_y_velocidad(Particulas &p, int posA){
    accX = p.fuerzasX[posA]/p.masas[posA];
    accY = p.fuerzasY[posA]/p.masas[posA];
    accZ = p.fuerzasZ[posA]/p.masas[posA];
    p.velocidadesX[posA]+= accX*time_step;
    p.velocidadesY[posA]+= accY*time_step;
    p.velocidadesZ[posA]+= accZ*time_step;
}

//x=x0+v*dt
void actualizar_posicion(Particulas &p, int posA){
    p.posicionesX[posA]+= p.velocidadesX[posA]*time_step;
    p.posicionesY[posA]+= p.velocidadesY[posA]*time_step;
    p.posicionesZ[posA]+= p.velocidadesZ[posA]*time_step;
    if (p.posicionesX[posA]<=0){
        p.posicionesX[posA]=0;
        p.velocidadesX[posA]=-p.velocidadesX[posA];
    }
    if (p.posicionesY[posA]<=0){
        p.posicionesY[posA]=0;
        p.velocidadesY[posA]=-p.velocidadesY[posA];
    }
    if (p.posicionesZ[posA]<=0){
        p.posicionesZ[posA]=0;
        p.velocidadesZ[posA]=-p.velocidadesZ[posA];
    }
    if (p.posicionesX[posA]>=size_enclosure){
        p.posicionesX[posA]=size_enclosure;
        p.velocidadesX[posA]=-p.velocidadesX[posA];
    }
    if (p.posicionesY[posA]>=size_enclosure){
        p.posicionesY[posA]=size_enclosure;
        p.velocidadesY[posA]=-p.velocidadesY[posA];
    }
    if (p.posicionesZ[posA]>=size_enclosure){
        p.posicionesZ[posA]=size_enclosure;
        p.velocidadesZ[posA]=-p.velocidadesZ[posA];
    }
    //reiniciar las fuerzas despues de actualizar la posicion
    p.fuerzasX[posA]=0;
    p.fuerzasY[posA]=0;
    p.fuerzasZ[posA]=0;

}
//crear el archivo final_config.txt
void final_config(Particulas &particulas){
    //crea el archivo en cmake-buiild-debug
    fstream file;
    file.open("final_config.txt",fstream::in | fstream::out | fstream::trunc);
    file<<fixed<<showpoint;
    file<<setprecision(3);
    file << size_enclosure<<" "<<time_step<<" "<<num_objects<<endl;
    for(int i = 0; i < num_objects; i++){

        file<<particulas.posicionesX[i]<<" "<<particulas.posicionesY[i]<<" "<<particulas.posicionesZ[i]<<" ";
        file<<particulas.velocidadesX[i]<<" "<<particulas.velocidadesY[i]<<" "<<particulas.velocidadesZ[i]<<" ";
        file<<particulas.masas[i]<<endl;

    }
    file.close();



}
int main(int argc, char* argv[]) {
    cout << "Has introducido " << argc-1 << " Parametros" << endl;
    //[----------------Control de Parámetros---------------------]
    if(argc != 6){
        cerr << "Número de parámetros Inválido" << endl;
        return -1;
    }
    ////dar valor a las variables globales
    num_objects = atoi(argv[1]);
    num_iteration = atoi(argv[2]);
    random_seed =atoi(argv[3]);
    size_enclosure = atof(argv[4]);
    time_step =atof(argv[5]);

    if(num_objects <= 0){
        cerr << "Número de partículas Inválido" << endl;
        cerr << "El número de partículas debe ser Mayor que 0" << endl;
        return -1;
    }else if(num_iteration <= 0){
        cerr << "Número de iteraciones de tiempo Inválido" << endl;
        cerr << "El número de iteraciones de tiempo debe ser Mayor que 0" << endl;
        return -1;
    }else if(random_seed <= 0){
        cerr << "Seed Inválida" << endl;
        cerr << "Seed debe ser Mayor que 0" << endl;
        return -1;
    }else if( size_enclosure<= 0){
        cerr << "Tamaño del cubo Inválido" << endl;
        cerr << "El tamaño del cubo debe ser Mayor que 0" << endl;
        return -1;
    }else if(time_step <= 0){
        cerr << "Intervalo de Tiempo Inválido" << endl;
        cerr << "El Intervalo de Tiempo debe ser Mayor que 0" << endl;
        return -1;
    }
    cout<<"Argumentos:\n"<<"num_objects: "<<num_objects<<"\nnum_iterations: "<<num_iteration<<"\nrandom_seed: "<<random_seed<<"\nsize_enclosure: "<<size_enclosure<<"\ntime_step: "<<time_step<<endl;
    //generacion de numeros aleatorios
    mt19937_64 generator(random_seed);
    uniform_real_distribution<double> dis{0.0, size_enclosure};
    normal_distribution<double> d{pow(10.0,21.0), pow(10.0, 15.0)};
    //crear la lista de particulas, darles valores iniciales y crear el init_config.txt
    Particulas particulas;
    generar_particulas(particulas,generator,dis,d);
    init_config(particulas);
    //por cada iteracion y particula, realizar las funciones anteriormente descritas
    for (int i=0;i<num_iteration;i++){
        for(int j=0;j<num_objects;j++){
            fuerza_gravitatoria(particulas,j);
            aceleracion_y_velocidad(particulas,j);
            actualizar_posicion(particulas,j);
            int pos = comprobar_colision(particulas, j);
            if(pos != -1){
                colision_particulas(particulas,j, pos);
            }

        }
    }
    //escribir el final_config.txt
    final_config(particulas);
    return 0;
}
