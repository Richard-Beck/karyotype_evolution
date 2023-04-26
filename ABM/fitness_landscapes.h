// generic fitness landscape that should be able to hold gaussian/polyh/random_field data
// and access the appropriate methods.
class fitness_landscape{
    public:
        //used in all landscapes
        list<vector<int>> peaks;
        fitness_landscape(list<vector<int>>& p){
            peaks = p;
        }
        virtual float get_fitness(vector<int>&)=0;

        //overloaded initializer signatures for various landscapes
        virtual void init(list<vector<int>>&,vector<float>&,vector<float>&){};
        virtual void init(vector<vector<float>>&,vector<float>&,vector<float>&){};
        virtual void init(list<vector<int>>&,vector<float>&, vector<float>&,vector<float>&, int){};
        virtual void init(list<vector<int>>&, float, float, float){};
        virtual void init(int, string){};
};

class krig_landscape: public fitness_landscape{
    public:
        vector<vector<float>> knots;
        vector<float> c;
        vector<float> d;
        // the p doesn't do anything for this landscape (although it could). The whole fitness landscape landscape needs reworked
        krig_landscape(list<vector<int>>& p):fitness_landscape(p) {}
        void init(vector<vector<float>>& k, vector<float>& cc,vector<float>& dd){
            knots=k;
            c=cc;
            d=dd;
        };
        float get_fitness(vector<int>&);

};

float krig_landscape::get_fitness(vector<int>& cn){
    float f=0, xx1=0;
    vector<float> xx0;// = {1.0};
    for(const auto& i:cn) xx0.push_back((float)i);

//    for(int i = 0; i<xx0.size(); i++){
  //      xx1+=d[i]*xx0[i];
    //}

     for(int i = 0; i<knots.size(); i++){
        float Di = 0;
        for(int j = 0; j<knots[i].size(); j++){
            Di+=pow((knots[i][j]-xx0[j]),2);
        }
        f+=c[i]*exp(-sqrt(Di));
     }
     f+=d[0];//xx1;
    return(f);
}

//Broken??
//float krig_landscape::get_fitness(vector<int>& cn){
    //float f=0, xx1=0;
    //vector<float> xx0 = {1.0};
    //for(const auto& i:cn) xx0.push_back((float)i);

    //for(int i = 0; i<xx0.size(); i++){
  //      xx1+=d[i]*xx0[i];
//    }

     //for(int i = 0; i<knots.size(); i++){
   //     float Di = 0;
 //       for(int j = 0; j<knots[i].size(); j++){
          //  Di+=pow((knots[i][j]-xx0[j+1]),2);
        //}
       // f+=c[i]*exp(-sqrt(Di));
     //}
    // f+=xx1;
  //  return(f);
//}

class gaussian_landscape: public fitness_landscape{
    public:
        vector<float> heights;
        vector<float> sigmas;
        gaussian_landscape(list<vector<int>>& p):fitness_landscape(p) {}
        void init(list<vector<int>>& p, vector<float>& h,vector<float>& s){
            peaks=p;
            heights=h;
            sigmas=s;
        };
        float get_fitness(vector<int>&);

};

float gaussian_landscape::get_fitness(vector<int>& cn){
    int i = 0;
    float ff = 0.0;
    for(const auto pk:peaks){
        float ffi = 0.0;
        for(int j =0; j<cn.size(); j++){
            ffi+=pow((cn[j]-pk[j])/sigmas[i],2);
        }
        ffi = heights[i]*exp(-ffi);
        ff = max(ff,ffi);
        i++;
    }
    return(ff);
}

class polyharmonic_landscape: public fitness_landscape{
    public:
        int k;
        vector<float> f;
        vector<float> w;
        vector<float> v;
        polyharmonic_landscape(list<vector<int>>& p):fitness_landscape(p) {

        }
        void init(list<vector<int>>& p, vector<float>& ff, vector<float>& ww,
                               vector<float>& vv, int kk){
            peaks=p;
            k=kk;
            f=ff;
            v=vv;
            w=ww;
        }
        float rbf(float r);
        float get_fitness(vector<int>&);

};

float polyharmonic_landscape::rbf(float r){
    float val;
    if(r>=1){
        if(k%2==0){
            val=pow(r,k)*log(r);
        }else{
            val = pow(r,k);
        }
    }else{
        if(k%2==0){
            val=pow(r,k-1)*log(pow(r,r));
        }else{
            val = pow(r,k);
        }
    }
    return val;
}

float polyharmonic_landscape::get_fitness(vector<int>& pt){
    vector<int> p2;
    p2.push_back(1);
    for(int i=0; i<pt.size(); i++) p2.push_back(pt[i]);

    float min_dist = 100.0;

    float fitness = 0;
    vector<float> r; // consider pre-allocating
    for(const auto knot_i:peaks){
        float d_i = 0;
        for(int j = 0; j<pt.size(); j++){
            d_i+=pow(knot_i[j]-pt[j],2);
        }
        min_dist=min(min_dist,sqrt((float)d_i));
        r.push_back(rbf(sqrt((float)d_i)));
    }
    for(int i = 0; i<r.size(); i++){
        fitness+=w[i]*r[i];
    }
    for(int i = 0; i<p2.size(); i++){
        fitness+=v[i]*p2[i];
    }
    fitness*=(1.0-pow(min_dist,10.0)/(1.0+pow(min_dist,10.0)));
    return fitness;
}

class random_landscape: public fitness_landscape{
    public:
        float wavelength;
        float scale;
        float centre;
        random_landscape(list<vector<int>>& p):fitness_landscape(p) {}
        void init(list<vector<int>>& p, float w, float s, float c){
            peaks = p;
            wavelength = w;
            scale = s;
            centre = c;
        };
        float get_fitness(vector<int>&);
};

float random_landscape::get_fitness(vector<int>&cn){
    float fitness = 0;
    for(const auto p:peaks){
        float d=0;
        for(int i =0; i<cn.size(); i++){
            d+=pow((float)cn[i]-p[i],2);
        }
        d=sqrt(d);
        fitness+=sin(d/wavelength);
    }

    fitness*=scale;
    fitness+=centre;
    return fitness;
}


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
    vector<vector<float>> pca_rot;
    vector<float> pca_cen;
    vector<float> t;
    int nchrom = 0;
    int nvar=0;
    int nterms=0;
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
    vector<float> y;

    for(int i = 0; i<nvar; i++){
        float yi=0.;
        for(int j = 0; j<z.size(); j++){
            float zj= z[j];
            yi+=(zj-pca_cen[j])*pca_rot[j][i];
        }
        //cout << yi << " ";
        y.push_back(yi);
    }
    //cout << endl;

    for(int i=0;i<y.size(); i++){
        y[i]= (y[i]-xc[i])/xs[i];
        //cout << y[i] << " ";
    }
    //cout << endl;
    int nt = 0;
    int d = y.size();
    if(nterms>1){
        for(int j = 0; j<d; j++){
            nt = j+1;
            wptr[j] = nt;
            ptab[nt][j] = ptab[nt][j]+1;
            t[nt] = y[j];
           // cout << t[nt] << endl;
        }
        //cout << "..." << endl;
        if(m>2){
            for(int k = 1; k<(m-1); k++){
                for(int j = 0; j<d; j++){

                    int bptr = wptr[j];
                    wptr[j] = nt+1;
                    int eptr = wptr[0]-1;
                    //cout << "bptr: " << bptr << " eptr: " << eptr << endl;
                    for(int tt=bptr; tt<(eptr+1);tt++){
                        nt++;
                        for(int jj = 0; jj<d; jj++){
                            ptab[nt][jj] = ptab[tt][jj];
                        }
                        ptab[nt][j]++;
                        t[nt]=y[j]*t[tt];
                        //cout << "tt: " <<tt << " t[tt]: "<<t[tt] <<" y[j]: "<< y[j] << endl;
                    }
                }
            }
        }
    }
    //cout << endl << t.size() << endl;
    //for(const auto& tit:t) cout << tit << " " << endl;

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
        if(words[0]=="cpca"){
            for(int i=1; i<words.size(); i++){
                pca_cen.push_back(stof(words[i]));
            }
        };
        if(words[0]=="knots"){
            int i = 1;
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
        if(words[0]=="rpca"){
            int i = 1;
            while(i<words.size()){
                int j = 0;
                vector<float> ri;
                while(j<nchrom){
                    j++;
                    ri.push_back(stof(words[i]));
                    i++;
                }
                pca_rot.push_back(ri);
            }
        };
    }
    nterms = binomialCoefficients((m+nvar-1),nvar);
    vector<int> wptr_v(nvar*m,0);
    vector<vector<int>> ptab_v(nterms,vector<int>(nvar, 0));
    vector<float> t_v(nterms,0);
    wptr=wptr_v;
    ptab = ptab_v;
    t=t_v;

}

struct hoc_landscape{

    float sigma =1;
    float mean =1.;
    float f0=1.0;
    map<vector<int>,float> visited;

    hoc_landscape() = default;
    hoc_landscape(float s, float m){
        sigma=s;
        mean=m;
        }
    float get_fitness(vector<int>&);
    float get_fitness(vector<int>&, float, mt19937&);
};

float hoc_landscape::get_fitness(vector<int>& cn){
    return f0;
}


float hoc_landscape::get_fitness(vector<int>& cn, float f, mt19937& gen){
    if (auto search = visited.find(cn); search != visited.end()){
            return(search->second);
            }else{ //if this is a new clone then generate a new entry in the map
                normal_distribution<> d{f+mean,sigma};
                f=d(gen);
                visited[cn]=f;
            }
    return f;
}



