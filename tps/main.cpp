#include <iostream>
#include <vector>
#include <float.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <list>
using namespace std;
float PI = 3.141592;

class fitness_landscape{
    public:
        //used in all landscapes
        list<vector<int>> peaks;
        fitness_landscape(list<vector<int>>& p){
            peaks = p;
        }
        virtual float get_fitness(vector<int>&)=0;

        //overloaded initializer signatures for various landscapes
        //virtual void init(list<vector<int>>&,vector<float>&,vector<float>&){};
        //virtual void init(list<vector<int>>&,vector<float>&, vector<float>&,vector<float>&, int){};
        //virtual void init(list<vector<int>>&, float, float, float){};
        virtual void init(int, string){};
};

int binomialCoefficients(int n, int k) {
   if (k == 0 || k == n)
   return 1;
   return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
}

float gamma_local(float x) {
  if (x < 0.) {
    float temp = 1.;
    while (x < 0.) {
      temp *= x;
      x++;
    }
    return tgamma(x)/temp;
  }
    return tgamma(x);
}

float rbs_const(float m,int d){
    float d2=d;
    d2/=2.;
    float c;
    if (d%2 == 0) {
            c= (pow(-1.,1.+m+d2) * pow(2.,1.-2*m) * pow(PI,-d2))/(tgamma(m)*gamma_local(m-d2+1.0));
        }else {
            c= (gamma_local(d2-m)*pow(2.,-2.*m)* pow(PI,-d2))/tgamma(m);
        }
    return c;
}

float rbf(float d, vector<float> par){
  d <- max(d,FLT_MIN);
  if(par[1]==0.){
    return pow(d,par[0]);
  }
  return log(d)*pow(d,par[0])/2;

}

float mrb(vector<float>& x1, vector<vector<float>>& x2,
                      vector<float>& C, vector<float>& par){

  float s0=0;
  for(int i = 0; i<x2.size(); i++){
    float d = 0.;
    for(int j = 0; j<x1.size();j++){
        d+=pow(x1[j]-x2[i][j],2.);
    }
    s0+=rbf(d,par)*C[i];
  }
  return s0;

}


class tps_landscape: public fitness_landscape{
    public:
    vector<float> xc;
    vector<float> xs;
    vector<float> dd;
    vector<float> c;
    vector<int> wptr;
    vector<vector<int>> ptab;
    vector<float> t;
    int nchrom = 0;
    int nterms=0;
    int nvar=0;
    vector<vector<float>> knots;
    int m;
    float p;
    tps_landscape(list<vector<int>>& pp):fitness_landscape(pp){}
    void init(int, string);
    float get_fitness(vector<int>&);
};

float tps_landscape::get_fitness(vector<int>& z){

    float f = 0.;
    t[0] = 1;
    vector<float> y(z.size(),0.);

    for(int i=0;i<z.size(); i++){
        float zi = z[i];
        y[i]= (zi-xc[i])/xs[i];
    }
    int nt = 0;
    int d = y.size();
    if(nterms>1){
        for(int j = 0; j<d; j++){
            nt = j+1;
            wptr[j] = nt;
            ptab[nt][j] = ptab[nt][j]+1;
            t[nt] = y[j];
        }
        if(m>2){
            for(int k = 1; k<(m-1); k++){
                    cout << k <<" ";
                for(int j = 0; j<d; j++){
                    int bptr = wptr[j];
                    wptr[j] = nt+1;
                    int eptr = wptr[0]-1;
                    for(int tt=(bptr-1); tt<eptr;tt++){
                        nt++;
                        for(int jj = 0; jj<d; jj++){
                            ptab[nt][jj] = ptab[tt][jj];
                        }
                        ptab[nt][j]++;
                        t[nt]=y[j]*t[tt];
                    }
                }
            }
        }
    }

    for(int i = 0; i<dd.size(); i++){
        f+=dd[i]*t[i];
    }

    float p2 = 0.;
    if(y.size()%2==0) p2 = 1.;
    vector<float> par = {p/2,p2};
    float m2 = ((float)y.size()+p)/2.;
    f+=mrb(y,knots,c,par)*rbs_const(m2,y.size());
    return(f);
};

void tps_landscape::init(int n, string path){
    nchrom=n;
    fstream fin;
    fin.open(path, ios::in);
    std::string tmp, row;
    vector<string> words;
    char delim = ',';

    while(std::getline(fin, row)){
        words.clear();
        stringstream tmp2(row);
        while (std::getline(tmp2, tmp, delim)) {
            words.push_back(tmp);
        }
        if(words[0]=="nvar") nvar=stoi(words[1]);
        if(words[0]=="m") m=stoi(words[1]);
        if(words[0]=="p") p=stof(words[1]);
        if(words[0]=="d"){
            for(int i=1; i<words.size(); i++){
                dd.push_back(stof(words[i]));
            }
        };
        if(words[0]=="c"){
            for(int i=1; i<words.size(); i++){
                c.push_back(stof(words[i]));
            }
        };
        if(words[0]=="xc"){
            for(int i=1; i<words.size(); i++){
                xc.push_back(stof(words[i]));
            }
        };
        if(words[0]=="xs"){
            for(int i=1; i<words.size(); i++){
                xs.push_back(stof(words[i]));
            }
        };
        if(words[0]=="knots"){
            int i = 1;
            cout << "nvar:" << nvar << endl;
            while(i<words.size()){
                int j = 0;
                vector<float> ki;
                while(j<nvar){
                    j++;
                    ki.push_back(stof(words[i]));
                    i++;
                }
                knots.push_back(ki);
            }
        };
    }
    nterms = binomialCoefficients((m+nchrom-1),nchrom);
    vector<int> wptr_v(nchrom*m,0);
    vector<vector<int>> ptab_v(nterms,vector<int>(nchrom, 0));
    vector<float> t_v(nterms,0);
    wptr=wptr_v;
    ptab = ptab_v;
    t=t_v;

}

int main()
{

    list<vector<int>> p;
    vector<int> z = {7,8};
    string path = "C:/Users/4473331/Documents/projects/008_birthrateLandscape/karyotype_evolution/ABM/landscapes/tst_tps.txt";
    tps_landscape l(p);
    l.init(2,path);
    cout <<"f= ";
    cout << l.get_fitness(z);

    return 0;
}
