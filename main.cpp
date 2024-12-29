#include<iostream>
#include<fstream>
#include<numeric>
#include <random>
#include <vector>
#include <queue>
#include <cstdint>
#include <cmath>
#include <array>
#include <stack>


using namespace std;

const int Np = 4 * 50 * 50 * 3;
int Lbox=50; // width/2

double pos[Np][3];
double forces[Np][3];
//uint16_t nearest[Np][3];
vector<uint16_t> nearest[Np];

// when calculating forces for particles, only run calculation if i_this < i_connected

int main(/*int argc=0, char** argv=nullptr*/){
    setQDs();
    srand((unsigned) time(NULL));


    apd=vector<vector<uint32_t>>(time_resolution,vector<uint32_t>(energy_resolution));
    pos p(0,0,0);
    int w_6=width/6;
    int w_23=width*2/3;
    for(uint32_t i=0; i<2000000; i++){
        p=pos((uint16_t)(w_6+rand()%w_23),(uint16_t)(w_6+rand()%w_23),(uint16_t)(rand()%layers));
        sim_particle(p,0);
    }
    save_as_csv(apd, "output.csv");
    return 0;
}