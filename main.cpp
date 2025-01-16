#include<iostream>
#include<fstream>
#include<numeric>
#include <random>
#include <vector>
#include <queue>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <array>
#include <stack>


using namespace std;
int Np;
int Lbox; // width

int ncell;
uint8_t num_cells;
int cell_w;
uint8_t max_z_cell;
vector<vector<vector<vector<uint32_t>>>> cells;

uint16_t d_size;
uint32_t *distances;
double d_res;

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
    double distance2(const vec3& p1) const { return (x - p1.x) * (x - p1.x) + (y - p1.y) * (y - p1.y) + (z - p1.z) * (z - p1.z);};
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

struct pos_r{
    pos_r(){i=0;r=0;}
    pos_r(uint32_t i1,double r1){i=i1;r=r1;}
    bool operator<=(const pos_r& p1) const{ return r<=p1.r;};
    bool operator>=(const pos_r& p1) const{ return r>=p1.r;};
    bool operator<(const pos_r& p1) const{ return r<p1.r;};
    bool operator>(const pos_r& p1) const{ return r>p1.r;};
    bool operator==(const pos_r& p1) const{ return r==p1.r;};
    uint32_t i; double r;
};

void remove_from_cell(uint8_t x,uint8_t y,uint8_t z,uint32_t i){
    auto it = std::lower_bound(cells[x][y][z].begin(), cells[x][y][z].end(), i);

    // Check if X is present
    if (it != cells[x][y][z].end() && *it == i) {
        cells[x][y][z].erase(it); // Remove the element
    }
}
void add_to_cell(uint8_t x,uint8_t y,uint8_t z,uint32_t i){
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
double LJ_sig = 2;
double max_force=2000;
bool too_close=false;
vec3 get_force(vec3 a, vec3 b, bool offset=false, double xa_off=0, double ya_off=0){
    vec3 r_v = b - a;
    if(offset){
        r_v.x-=xa_off;
        r_v.y-=ya_off;
        //cout<<"offset = "<<xa_off<<", "<<ya_off<<"\n";
    }
    double r_mag = sqrt(r_v.x*r_v.x+r_v.y*r_v.y+r_v.z*r_v.z);
    if(r_mag==0){
        a.print();cout<<" is for some reason applying force to itself\n";
        too_close= true;
        return {0,0,0};
    }
    double f = - 24 * LJ_eps * (pow((LJ_sig / r_mag),7)) * (1 - 2 * pow(LJ_sig / r_mag,6)) / (LJ_sig);
    if(r_mag<1.5){
        a.print();cout<<" and ";b.print();cout<<" have r = "<<r_mag<<", f = "<<f<<"\n";
        too_close= true;
        if(f>100){f=100;}else if(f<-100){ f=-100;}
    }
    r_v.scale(f/(r_mag*r_mag));
    return r_v;

}

void force_boundary(uint8_t x0,uint8_t y0,uint8_t z0,uint8_t x1,uint8_t y1,uint8_t z1, vec3 forces[], vec3 pos[], cell_pos cpos[]){
    double x0_offset=0;
    double y0_offset=0;
    if(x0==0 && x1+1==num_cells){
        x0_offset=2*Lbox;
    }else if(x0+1==num_cells && x1==0){
        x0_offset=-2*Lbox;
    }
    if(y0==0 && y1+1==num_cells){
        y0_offset=2*Lbox;
    }else if(y0+1==num_cells && y1==0){
        y0_offset=-2*Lbox;
    }
    for (auto & i : cells[x0][y0][z0]){
        for (auto & j : cells[x1][y1][z1]){
            vec3 force= get_force(pos[i], pos[j], true,x0_offset,y0_offset);
            forces[i].subtract(force);
            forces[j].add(force);
            if(too_close){
                cout<<i<<"-"<<j<<" ";cpos[i].print();cout<<" - ";cpos[j].print();cout<<" \n";
                cout<<i<<"-"<<j<<" ("<<unsigned(x0)<<","<<unsigned(y0)<<","<<unsigned(z0)<<") - ("<<unsigned(x1)<<","<<unsigned(y1)<<","<<unsigned(z1)<<"), boundary \n";
                too_close= false;
            }
        }
    }

}

void force_cells(uint8_t x0,uint8_t y0,uint8_t z0,uint8_t x1,uint8_t y1,uint8_t z1, vec3 forces[], vec3 pos[], cell_pos cpos[]){
    if(x1>num_cells-1||y1>num_cells-1){
        if(x1>num_cells-1){
            x1=(x1==num_cells)?0:num_cells-1;
        }
        if(y1>num_cells-1){
            y1=(y1==num_cells)?0:num_cells-1;
        }
        //cout<<"boundary from "<<unsigned(x0)<<","<<unsigned(y0)<<","<<unsigned(z0)<<" to "<<unsigned(x1)<<","<<unsigned(y1)<<","<<unsigned(z0)<<"\n";
        force_boundary(x0,y0,z0,x1,y1,z1,forces,pos,cpos);
        return;
    }
    for (auto & i : cells[x0][y0][z0]){
        for (auto & j : cells[x1][y1][z1]){
            vec3 force= get_force(pos[i], pos[j]);
            forces[i].subtract(force);
            forces[j].add(force);
            if(too_close){
                cout<<i<<"-"<<j<<" (";cpos[i].print();cout<<") - (";cpos[j].print();cout<<") \n";
                cout<<i<<"-"<<j<<" ("<<unsigned(x0)<<","<<unsigned(y0)<<","<<unsigned(z0)<<") - ("<<unsigned(x1)<<","<<unsigned(y1)<<","<<unsigned(z1)<<") \n";
                too_close= false;
            }
        }
    }
}

void force_self(uint8_t x,uint8_t y,uint8_t z, vec3 forces[], vec3 pos[], cell_pos cpos[]) {
    uint8_t sz=uint8_t(cells[x][y][z].size());
    for (uint8_t i=0; i<sz; i++) {
        uint32_t i_i=cells[x][y][z][i];
        for(uint8_t j=i+1;j<sz;j++){
            uint32_t j_i=cells[x][y][z][j];
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


void get_nns(uint32_t i, const vec3 pos[], const cell_pos pos_c[], ofstream *file){
    uint8_t cx0=pos_c[i].x;
    uint8_t cy0=pos_c[i].y;
    uint8_t cz0=pos_c[i].z;
    uint8_t cx_max=(cx0+1==ncell)?cx0:cx0+1;
    uint8_t cy_max=(cy0+1==ncell)?cy0:cy0+1;
    uint8_t cz_max=(cz0+1==ncell)?cz0:cz0+1;
    vector<pos_r> closest;
    for(uint8_t cx=(0<cx0)?cx0-1:cx0; cx<=cx_max; cx++){
        for(uint8_t cy=(0<cy0)?cy0-1:cy0; cy<=cy_max; cy++) {
            for(uint8_t cz=(0<cz0)?cz0-1:cz0; cz<=cz_max; cz++) {
                for (uint32_t j:cells[cx][cy][cz]){
                    if(j==i){continue;}
                    double r_mag = pos[i].distance2(pos[j]);
                    if(r_mag<16){
                        uint16_t r_in=(uint16_t)(sqrt(r_mag)/d_res);
                        //cout<<r_in<<"/"<<d_size<<" ";
                        if(r_in<d_size){
                            distances[r_in]++;
                        }
                    }
                    pos_r pr=pos_r(j,r_mag);

                    auto it = std::lower_bound(closest.begin(),closest.end(), pr);

                    closest.insert(it, pr);
                    //cout<<pr.i<<" ("<<pr.r<<"|"<<r_mag<<"), ";

                }
            }
        }
    }
    //cout<<i<<": ";
    //cout<<"["<<closest.size()<<"] ";
    int num=0;
    for(uint32_t k=0; k<min(closest.size(),size_t(12));k++){
        //if(closest[k].r>8){break;}
        //if(k<4)
        //    cout<<closest[k].i<<" ("<<closest[k].r<<"), ";

        *file<<","<<closest[k].i;
        num++;
    }
    //cout<<"\n";
    if(num!=12){
        cout<<"\nERROR, not 12 NNS, instead there are "<<num<<"\n";
    }
    //cout<<num<<" ";

    //cout<<"\n";
    *file<<"\n";
}


void save_as_csv(const vec3 pos[], const cell_pos pos_c[], const std::string& filename) {
    auto start = chrono::high_resolution_clock::now();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    for (int n=0;n<Np;n++) {
        file<<pos[n].x<<","<<pos[n].y<<","<<pos[n].z;//<<"\n";
        get_nns(n,pos,pos_c,&file);
    }



    auto end = chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    file.close();
    std::cout << "CSV saved to " << filename <<" in "<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<" ms\n";

}

void save_dist_plot(const std::string& filename) {

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
        return;
    }
    for (uint16_t d=0;d<d_size;d++) {
        file<<(d_res*d)<<","<<distances[d]<<"\n";
    }

    file.close();

}


int main(/*int argc=0, char** argv=nullptr*/){
    srand((unsigned) time(NULL));

    d_res=.003;
    Lbox=210;
    cell_w=3;
    ncell=(2*Lbox)/cell_w;
    max_z_cell=(uint8_t)ncell-1;
    cells = vector<vector<vector<vector<uint32_t>>>>(ncell,vector<vector<vector<uint32_t>>>(ncell,vector<vector<uint32_t>>(ncell,vector<uint32_t>())));
    Np = (int)(Lbox * Lbox * 3.2);
    d_size=(uint16_t)ceil(4/d_res);
    distances=new uint32_t[d_size];
    fill_n(distances,d_size,0);
    //cout<<"distances low: "<<distances[0]<<", "<<distances[1]<<", "<<distances[2]<<", "<<distances[3]<<"\n";

    vec3 *pos = new vec3[Np];
    cell_pos *pos_c = new cell_pos[Np];
    vec3 *vel = new vec3[Np];
    vec3 *forces = new vec3[Np];
    set_init(pos,pos_c);
    //double total_dx=100;
    int iter=0;
    double dt=.05;
    double avg_E=1000;
    num_cells=(uint8_t)ncell;
    auto start = chrono::high_resolution_clock::now();
    double last_E=1;
    double last_az=1;
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
        uint8_t highest_z=0;
        for(uint8_t z0=0; z0<=max_z_cell; z0++){
            for(uint8_t x0=0;x0<num_cells;x0++){
                //uint8_t x1_max=(x0==num_cells-1)?x0:x0+1;
                for(uint8_t y0=0; y0<num_cells; y0++){
                    //uint8_t y1_max=(y0==num_cells-1)?y0:y0+1;
                    if(cells[x0][y0][z0].empty()){continue;}
                    if(z0>highest_z){
                        highest_z=z0;
                    }
                    if(cells[x0][y0][z0].size()>1){
                        force_self(x0,y0,z0,forces,pos,pos_c);
                    }
                    force_cells(x0,y0,z0,x0+1,y0,z0,forces,pos,pos_c);
                    uint8_t x1=x0-1;
                    for(int x1o=-1; x1o<=1; x1o++) {
                        force_cells(x0,y0,z0,x1,y0+1,z0,forces,pos,pos_c);
                        x1++;
                    }
                    if(z0+1<=max_z_cell){
                        uint8_t y1=y0-1;
                        for(int y1o=-1; y1o<=1; y1o++) {
                            x1=x0-1;
                            for(int x1o=0; x1o<=2; x1o++) {
                                force_cells(x0,y0,z0,x1,y1,z0+1,forces,pos,pos_c);
                                x1++;
                            }
                            y1++;
                            /*for(uint8_t x1o=(x0==0)?0:x0-1; x1<=y0+1; x1++) {
                                force_cells(x0,y0,z0,x1,y1,z0+1,forces,pos,pos_c);
                            }*/
                        }
                    }
                    /*if(x0+1<num_cells){
                        force_cells(x0,y0,z0,x0+1,y0,z0,forces,pos,pos_c);
                    } else{
                        force_boundary(num_cells-1,y0,z0,0,y0,z0,forces,pos,pos_c);
                    }
                    if(y0+1<num_cells){
                        for(uint8_t x1=(x0==0)?0:x0-1; x1<=x1_max; x1++) {
                            force_cells(x0,y0,z0,x1,y0+1,z0,forces,pos,pos_c);
                        }
                    } else{
                        if(x0==0){ // x0 at left border
                            force_boundary(0,num_cells-1,z0,num_cells-1,0,z0,forces,pos,pos_c);
                        }else if(x0+1==num_cells){//x0 at right border
                            force_boundary(num_cells-1,num_cells-1,z0,0,0,z0,forces,pos,pos_c);

                        }
                        for(uint8_t x1=(x0==0)?0:x0-1; x1<=x1_max; x1++) {
                            force_boundary(x0,num_cells-1,z0,x1,0,z0,forces,pos,pos_c);
                        }
                    }
                    if(z0+1<=max_z_cell){
                        for(uint8_t y1=(y0==0)?0:y0-1; y1<=y1_max; y1++) {
                            for(uint8_t x1=(x0==0)?0:x0-1; x1<=x1_max; x1++) {
                                force_cells(x0,y0,z0,x1,y1,z0+1,forces,pos,pos_c);
                            }
                        }
                    }*/
                    //highest_z=z0;
                }
            }

        }
        max_z_cell=(uint8_t)highest_z;
        double avg_z=0;
        //double E_i=0;
        for(int i=0; i<Np;i++){
            forces[i].scale(dt);
            vel[i].add(forces[i]);
            pos[i].add_scaled(vel[i],dt);
            //E_i=vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z;
            avg_E+=vel[i].x*vel[i].x+vel[i].y*vel[i].y+vel[i].z*vel[i].z;
            /*if(sqrt(E_i)*dt>.5){
                cout<<i<<" is moving with vel "<<sqrt(E_i)*dt<<" - TOO FAST\n";
            }*/
            if(abs(pos[i].x) >= Lbox) {
                if(pos[i].x>Lbox){
                    pos[i].x-=Lbox*2;
                } else{
                    pos[i].x+=Lbox*2;

                }
            } if(abs(pos[i].y) >= Lbox ) {
                if(pos[i].y>Lbox){
                    pos[i].y-=Lbox*2;
                } else{
                    pos[i].y+=Lbox*2;
                }
            } if(abs(pos[i].z) >= Lbox - 1) {
                pos[i].z-=vel[i].z*dt;
                vel[i].z = -.8 * vel[i].z;
            }
            pos_c[i].update(pos[i],i);
            vel[i].scale(.95);
            forces[i].resetZ(1);
            avg_z+=pos[i].z;
        }
        if(iter%20==0){
            auto end = chrono::high_resolution_clock::now();
            avg_E=avg_E/Np;; avg_z=avg_z/Np+Lbox;
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            cout<<iter<<": avg z = "<<avg_z<<", E = "<<avg_E;
            if(iter>50){
                cout<<" ("<<round(100*(avg_E/last_E))<<"%)";
                if (avg_E/last_E>.98 && avg_E/last_E<1.1 && last_az-avg_z<.01){
                    dt=dt/2;
                    if(avg_E<.0001){
                        iter=5000;
                    }
                }
            }
            cout<<", highest z = "<<unsigned(max_z_cell)<<", dt = "<<dt;
            cout<<", ("<<std::chrono::duration_cast<std::chrono::milliseconds>(duration).count()<<" ms)\n";
            start=end; last_E=avg_E; last_az=avg_z;
        }
        iter++;
    }
    /*for(uint8_t z0=0; z0<=max_z_cell; z0++){
        for(uint8_t x0=0;x0<num_cells;x0++){
            for(uint8_t y0=0; y0<num_cells; y0++){
                cout<<"("<<unsigned(x0)<<","<<unsigned(y0)<<","<<unsigned(z0)<<"): ";
                for(uint32_t i:cells[x0][y0][z0]){
                    cout<<i<<", ";
                }
                cout<<"\n";
            }
        }

    }*/
    //cout<<"distances low: "<<distances[0]<<", "<<distances[1]<<", "<<distances[2]<<", "<<distances[3]<<"\n";
    string name="dots_p5_N"+to_string(Np)+"_w"+to_string(Lbox)+".csv";
    save_as_csv(pos, pos_c,name);

    /*for(uint16_t d=0; d<d_size; d++){
        cout<<d<<"-"<<distances[d]<<" | ";
    }*/
    string dists="distances_pf_N"+to_string(Np)+"_w"+to_string(Lbox)+"_95_1g.csv";
    save_dist_plot(dists);
    /*for(int i=0; i<Np;i++){
        pos[i].print();
        cout<<"\n";
    }*/

    //vector<vec> vel;
    //double forces[Np][3];
    //double forces[Np][3];


    delete [] distances;
    delete [] pos;
    delete [] pos_c;
    delete [] vel;
    delete [] forces;
    return 0;
}