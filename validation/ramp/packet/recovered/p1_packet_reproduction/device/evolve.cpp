#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <omp.h>
using namespace std;
// Cartesian positive y,z quadrant, reflected in between the innermost layers.
const int N=120,H=28,L=3*N*H*H+N;
int idx(int i,int j,int k,int component){return N+3*((i*H+j)*H+k)+component;}
int main(int argc,char **argv){
 if(argc!=2)return 2; string out=argv[1];
 vector<map<int,double>> entries(L);
 for(int i=0;i<N-1;i++){
   entries[i][i]-=1.69; entries[i][i+1]+=1.69; entries[i+1][i]+=1.69; entries[i+1][i+1]-=1.69;
 }
 for(int i=0;i<N;i++)for(int j=0;j<H;j++)for(int k=0;k<H;k++){
   for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
     const int n=dx*dx+dy*dy+dz*dz; if(n!=1 && n!=2)continue;
     if(i+dx<0 || i+dx>=N || j+dy>=H || k+dz>=H)continue;
     int sx=1,sy=(j+dy<0)?-1:1,sz=(k+dz<0)?-1:1;
     int y=sy<0?0:j+dy,z=sz<0?0:k+dz;
     int dir[]={dx,dy,dz}, reflection[]={sx,sy,sz};
     for(int l=0;l<3;l++)if(dir[l])for(int m=0;m<3;m++)if(dir[m]){
       double w=double(dir[l]*dir[m])/n;
       entries[idx(i,j,k,l)][idx(i,j,k,m)]-=w;
       entries[idx(i,j,k,l)][idx(i+dx,y,z,m)]+=w*reflection[m];
     }
   }
 }
 for(int i=60;i<=65;i++)for(int e=0;e<=1;e++){
   int row=i+e, vec[]={1-2*e,1,1};
   double g=((i==60||i==65)?.1:.2)/3;
   entries[row][row]-=4*g;
   for(int l=0;l<3;l++){
      int b=idx(i,0,0,l);
      entries[row][b]+=4*g*vec[0]*vec[l];entries[b][row]+=g*vec[l]*vec[0];
      for(int m=0;m<3;m++)entries[b][idx(i,0,0,m)]-=g*vec[l]*vec[m];
   }
 }
 vector<int> col,p(L+1);vector<double> val;
 for(int i=0;i<L;i++){p[i]=col.size();for(auto item:entries[i])if(abs(item.second)>1e-14){col.push_back(item.first);val.push_back(item.second);}}
 p[L]=col.size();entries.clear();entries.shrink_to_fit();
 vector<double> x(L),v(L),a(L);for(int j=0;j<N;j++)if(!(cin>>x[j]))return 3;for(int j=0;j<N;j++)if(!(cin>>v[j]))return 4;
 ofstream snap(out+"/snap.bin",ios::binary),line(out+"/line.bin",ios::binary),corners(out+"/corner.bin",ios::binary);
 double dt=.04;
 auto forces=[&](){
  #pragma omp parallel for schedule(static)
  for(int i=0;i<L;i++){double sum=0;for(int j=p[i];j<p[i+1];j++)sum+=val[j]*x[col[j]];a[i]=sum;}
 };
 forces();
 for(int n=0;n<=1425;n++){
   line.write((char*)x.data(),N*sizeof(double));line.write((char*)v.data(),N*sizeof(double));
   for(int iq=60;iq<=65;iq++) {corners.write((char*)(x.data()+idx(iq,0,0,0)),3*sizeof(double));corners.write((char*)(v.data()+idx(iq,0,0,0)),3*sizeof(double));}
   if(n%5==0 && n>=1200){
      snap.write((char*)x.data(),L*sizeof(double));snap.write((char*)v.data(),L*sizeof(double));
   }
   if(n==1425)break;
   #pragma omp parallel for schedule(static)
   for(int i=0;i<L;i++){v[i]+=.5*dt*a[i];x[i]+=dt*v[i];}
   forces();
   #pragma omp parallel for schedule(static)
   for(int i=0;i<L;i++)v[i]+=.5*dt*a[i];
 }
 return !(line && snap && corners);
}
