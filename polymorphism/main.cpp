#include <iostream>
#include <list>
#include <vector>
#include <math.h>

using namespace std;
#include "fitness_landscape.h"

class A{

    public:
        vector<int> v;
        A(vector<int>& a){
            v=a;
        }
        virtual void printv(){
            for(const auto & x:v) cout << " " << x;
            cout << endl;
        }
        virtual void reset(vector<int>& a, vector<int>& b)=0;
};

class B: public A{

    public:
        vector<int> w;
        B(vector<int>& a, vector<int>& b):A(a){
            w=b;
        }
        void reset(vector<int>& a, vector<int>& b){
            v=a;
            w=b;
        }
        void printv(){
            if(w.size()>v.size()){
                for(const auto & x:w) cout << " " << x;
            }else{
                for(const auto & x:v) cout << " " << x;
            }

            cout << endl;
        }
};

int main() {

    vector<int> v{2,2,2};
    vector<int> w{3,3,3,3};
    vector<int> w2{3,3,3,3,5};
    A *a;

    int c =0;
    cin >> c;
    B b(v,w);
    a = &b;
    if(c==1){
        a->reset(v,w2);

    }
    a->printv();

    list<vector<int>> pks;


    fitness_landscape *f;
    gaussian_landscape g(pks);


    vector<int> p1 = {1,1};
    vector<int> p2 = {2,2};


    pks.push_back(p2);
    pks.push_back(p1);

    vector<float> sigmas = {3.0,3.0};
    vector<float> heights = {1.0,1.0};

    if(c==1){
        f=&g;
        f->init(pks,sigmas,heights);

    }

    cout << f->get_fitness(p1);


    return 0;
}
