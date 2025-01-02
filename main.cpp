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
int Np;
int Lbox; // width

int ncell;
int cell_w;
uint8_t max_z_cell;
vector<vector<vector<vector<uint16_t>>>> cells;

struct vec3{
    vec3(){ x=0;y=0;z=0;}
    vec3(double x, double y, double z): x{x}, y{y}, z{z}{}
    vec3(vec3 &p0): x{p0.x}, y{p0.y}, z{p0.z}{}
    vec3(vec3 const &vec1) : x{vec1.x}, y{vec1.y}, z{vec1.z}{}
    void set_cell(){

    }
    void print(){cout << "(" << x << "," << y << "," << z << ") ";}
    void print_sub(vec3 &p1){cout << "(" << x - p1.x << "," << y - p1.y << "," << z - p1.z << ") ";}
    double get_dist(vec3 &p1){return pow(x - p1.x, 2) + pow(y - p1.y, 2) + pow(z - p1.z, 2);}
    void scale(double k){x=x*k;y=y*k;z=z*k;}
    void operator=(const vec3& p1){ x=p1.x; y=p1.y; z=p1.z;};
    void add(const vec3& p1){ x+=p1.x; y+=p1.y; z+=p1.z;};
    void add_scaled(const vec3& p1, double k){ x+= k * p1.x; y+= k * p1.y; z+= k * p1.z;};
    void subtract(const vec3& p1){ x-=p1.x; y-=p1.y; z-=p1.z;};
    void resetZ(double grav){x=0;y=0;z=-grav;}
    vec3 operator+(vec3 const& obj){
        vec3 res(x + obj.x, y + obj.y, z + obj.z);
        return res;
    }
    vec3 operator-(vec3 const& obj){
        vec3 res(x - obj.x, y - obj.y, z - obj.z);
        return res;
    }
    double x; double y; double z;
};

void remove_from_cell(uint8_t x,uint8_t y,uint8_t z,uint16_t i){
    auto it = std::lower_bound(cells[x][y][z].begin(), cells[x][y][z].end(), i);

    // Check if X is present
    if (it != cells[x][y][z].end() && *it == i) {
        cells[x][y][z].erase(it); // Remove the element
    }
}
void add_to_cell(uint8_t x,uint8_t y,uint8_t z,int i){
    auto it = std::lower_bound(cells[x][y][z].begin(), cells[x][y][z].end(), i);
    cells[x][y][z].insert(it, i);
    if(z>max_z_cell){
        max_z_cell=z;
    }
}

struct cell_pos{
    cell_pos(){x=0;y=0;z=0;};
    cell_pos(vec3 p,int i){
        x=(uint8_t)floor((p.x+Lbox)/cell_w);
        y=(uint8_t)floor((p.y+Lbox)/cell_w);
        z=(uint8_t)floor((p.z+Lbox)/cell_w);
        add_to_cell(x,y,z,i);
    }
    void set(vec3 p,int i){
        x=(uint8_t)floor((p.x+Lbox)/cell_w);
        y=(uint8_t)floor((p.y+Lbox)/cell_w);
        z=(uint8_t)floor((p.z+Lbox)/cell_w);
        add_to_cell(x,y,z,i);
    }
    void update(vec3 p,int i){
        uint8_t x1=(uint8_t)floor((p.x+Lbox)/cell_w);
        uint8_t y1=(uint8_t)floor((p.y+Lbox)/cell_w);
        uint8_t z1=(uint8_t)floor((p.z+Lbox)/cell_w);
        if(x1!=x||y1!=y||z1!=z){
            remove_from_cell(x,y,z,i);
            add_to_cell(x1,y1,z1,i);
            x=x1;y=y1;z=z1;
        }
    }
    void print(){cout << "(" << unsigned(x) << "," << unsigned(y) << "," << unsigned(z) << ") ";}
    uint8_t x; uint8_t y;uint8_t z;
};

//uint16_t nearest[Np][3];

// when calculating forces for particles, only run calculation if i_this < i_connected
void set_init(vec3 pos[],cell_pos pos_c[]){
    double x1= -Lbox;
    int y1= -Lbox;
    int z1= -Lbox;
    int step=3;
    for(int i=0; i<Np;i++){
        pos[i].x = x1+1+(step-2)*((rand()%1000)/1000.0);
        pos[i].y = y1+1+(step-2)*((rand()%1000)/1000.0);
        pos[i].z = z1+1+(step-2)*((rand()%1000)/1000.0);
        pos_c[i].set(pos[i],i);
        x1=pos[i].x+1;
        if(x1+step>=Lbox){
            x1=-Lbox;
            y1+=step;
            if(y1+step>=Lbox){
                y1=-Lbox;
                z1+=step;
                if(z1+step>=Lbox){
                    cout<<"can't fit all in box\n";
                    exit(1);
                }

            }
        }
        //cout<<"\n";
        //pos[i].print();
    }
}

double LJ_eps = .1;
double LJ_sig = 1.8;
double max_force=2000;
bool too_close=false;
vec3 get_force(vec3 a, vec3 b){
    vec3 r_v = b - a;
    double r_mag = sqrt(r_v.x*r_v.x+r_v.y*r_v.y+r_v.z*r_v.z);
    double f = - 24 * LJ_eps * (pow((LJ_sig / r_mag),8)) * (1 - 2 * pow(LJ_sig / r_mag,6)) / (LJ_sig*LJ_sig);
    if(r_mag<1.5){
        a.print();cout<<" and ";b.print();cout<<" have r = "<<r_mag<<", f = "<<f<<"\n";
        too_close= true;
        if(f>100){f=100;}else if(f<-100){ f=-100;}
    }
    r_v.scale(f/(r_mag*r_mag));
    /*if(r_v.x>max_force)
        r_v.x=max_force;
    else if(r_v.x<-max_force)
        r_v.x=-max_force;
    if(r_v.y>max_force)
        r_v.y=max_force;
    else if(r_v.y<-max_force)
        r_v.y=-max_force;
    if(r_v.z>max_force)
        r_v.z=max_force;
    else if(r_v.z<-max_force)
        r_v.z=-max_force;*/

    return r_v;

}

void force_cells(uint8_t x0,uint8_t y0,uint8_t z0,uint8_t x1,uint8_t y1,uint8_t z1, vec3 forces[], vec3 pos[], cell_pos cpos[]){
    for (auto & i : cells[x0][y0][z0]){
        for (auto & j : cells[x1][y1][z1]){
            vec3 force= get_force(pos[i], pos[j]);
            forces[i].subtract(force);
            forces[j].add(force);
            if(too_close){
                cout<<i<<"-"<<j<<" (";cpos[i].print();cout<<") - (";cpos[j].print();cout<<") \n";
                too_close= false;
            }
        }
    }
}

void force_self(uint8_t x,uint8_t y,uint8_t z, vec3 forces[], vec3 pos[], cell_pos cpos[]) {
    uint8_t sz=cells[x][y][z].size();
    for (uint8_t i=0; i<sz; i++) {
        uint16_t i_i=cells[x][y][z][i];
        for(uint8_t j=i+1;j<sz;j++){
            uint16_t j_i=cells[x][y][z][j];
            vec3 force= get_force(pos[i_i], pos[j_i]);
            forces[i_i].subtract(force);
            forces[j_i].add(force);

            if(too_close){
                cout<<i_i<<"-"<<j_i<<"(self) (";cpos[i_i].print();cout<<") - (";cpos[j_i].print();cout<<") \n";
                too_close= false;
            }

        }
    }
}

void save_as_csv(const vec3 pos[], const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    for (size_t n=0;n<Np;n++) {
        file<<pos[n].x<<","<<pos[n].y<<","<<pos[n].z<<"\n";
    }
    file.close();
    std::cout << "CSV saved to " << filename << std::endl;
}


int main(/*int argc=0, char** argv=nullptr*/){
    srand((unsigned) time(NULL));

    //ncell=10;
    Lbox=18;
    cell_w=3;
    ncell=(2*Lbox)/cell_w;
    max_z_cell=(uint8_t)ncell-1;
    cells = vector<vector<vector<vector<uint16_t>>>>(ncell,vector<vector<vector<uint16_t>>>(ncell,vector<vector<uint16_t>>(ncell,vector<uint16_t>())));
    Np = 20 * 20 * 4;

    vec3 *pos = new vec3[Np];
    cell_pos *pos_c = new cell_pos[Np];
    vec3 *vel = new vec3[Np];
    vec3 *forces = new vec3[Np];
    set_init(pos,pos_c);
    //double total_dx=100;
    int iter=0;
    double dt=.05;
    double avg_E=1000;
    uint8_t num_cells=(uint8_t)ncell;
    while(iter<5000){
        avg_E=0;
        /*for(int i=0; i<Np;i++){
            for(int j=i+1; j<Np;j++){
                vec3 force= get_force(pos[i], pos[j]);
                forces[i].subtract(force);
                forces[j].add(force);
            }
            //cout<<"particle ";pos[i].print();cout<<" has force ";forces[i].print();cout<<"\n";
        }*/
        //uint8_t highest_z=0;
        for(uint8_t z0=0; z0<=max_z_cell; z0++){
            for(uint8_t x0=0;x0<num_cells;x0++){
                uint8_t x1_max=(x0==num_cells-1)?x0:x0+1;
                for(uint8_t y0=0; y0<num_cells; y0++){
                    uint8_t y1_max=(y0==num_cells-1)?y0:y0+1;
                    if(cells[x0][y0][z0].empty()){continue;}
                    if(cells[x0][y0][z0].size()>1){
                        force_self(x0,y0,z0,forces,pos,pos_c);
                    }
                    if(x0+1<num_cells){
                        force_cells(x0,y0,z0,x0+1,y0,z0,forces,pos,pos_c);
                    }
                    if(y0+1<num_cells){
                        for(uint8_t x1=(x0==0)?0:x0-1; x1<=x1_max; x1++) {
                            force_cells(x0,y0,z0,x1,y0+1,z0,forces,pos,pos_c);
                        }
                    }
                    if(z0+1<=max_z_cell){
                        for(uint8_t y1=(y0==0)?0:y0-1; y1<=y1_max; y1++) {
                            for(uint8_t x1=(x0==0)?0:x0-1; x1<=x1_max; x1++) {
                                force_cells(x0,y0,z0,x1,y1,z0+1,forces,pos,pos_c);
                            }
                        }
                    }
                    //highest_z=z0;
                }
            }

        }
        //max_z_cell=(uint8_t)highest_z;
        double avg_z=0;
        for(int i=0; i<Np;i++){
            forces[i].scale(dt);
            vel[i].add(forces[i]);
            pos[i].add_scaled(vel[i],dt);
            avg_E+=vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z;
            if(abs(pos[i].x) >= Lbox - 1) {
                pos[i].x-=vel[i].x*dt;
                //pos[i].x = (pos[i].x<0) ? -Lbox+1 : Lbox-1;
                vel[i].x = -.8 * vel[i].x;
            } if(abs(pos[i].y) >= Lbox - 1) {
                pos[i].y-=vel[i].y*dt;
                //pos[i].y = (pos[i].y<0) ? -Lbox+1 : Lbox-1;
                vel[i].y = -.8 * vel[i].y;
            } if(abs(pos[i].z) >= Lbox - 1) {
                pos[i].z-=vel[i].z*dt;
                //pos[i].z = (pos[i].z<0) ? -Lbox+1 : Lbox-1;
                vel[i].z = -.8 * vel[i].z;
            }
            pos_c[i].update(pos[i],i);
            vel[i].scale(.6);
            forces[i].resetZ(2);
            avg_z+=pos[i].z;
        }
        if(iter%50==0)
            cout<<iter<<": avg z = "<<avg_z/Np<<", E = "<<avg_E/Np<<"\n";
        iter++;
    }

    string name="output2_N"+to_string(Np)+"_w"+to_string(Lbox)+".csv";
    save_as_csv(pos, name);
    /*for(int i=0; i<Np;i++){
        pos[i].print();
        cout<<"\n";
    }*/

    //vector<vec> vel;
    //double forces[Np][3];
    //double forces[Np][3];


    delete [] pos;
    delete [] pos_c;
    delete [] vel;
    delete [] forces;
    return 0;
}