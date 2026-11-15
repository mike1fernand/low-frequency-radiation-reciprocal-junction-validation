#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>

static inline size_t idx3(int i,int j,int k,int ny,int nz){return ((size_t)i*ny+j)*nz+k;}
struct Energies{double Eq=0,Eb=0,Ec=0;};
struct Shell{double union_total=0; std::array<double,6> face{}; double node_sum=0;};
const int NB=9;
const int B[NB][3]={{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,-1,0},{1,0,1},{1,0,-1},{0,1,1},{0,1,-1}};

void force_all(const std::vector<double>&q,const std::vector<double>&u,const std::vector<double>&g,double Kq,int nx,int ny,int nz,std::vector<double>&fq,std::vector<double>&fu){
  #pragma omp parallel for schedule(static)
  for(int i=0;i<nx;i++){double f=0; if(i>0)f+=Kq*(q[i-1]-q[i]); if(i+1<nx)f+=Kq*(q[i+1]-q[i]); fq[i]=f;}
  const size_t nc=(size_t)nx*ny*nz;
  #pragma omp parallel for schedule(static)
  for(long long pp=0;pp<(long long)nc;pp++){
    int k=pp%nz; long long rr=pp/nz; int j=rr%ny, i=rr/ny;
    double ux=u[3*pp],uy=u[3*pp+1],uz=u[3*pp+2],f0=0,f1=0,f2=0;
    for(int b=0;b<NB;b++){
      double L=std::sqrt(double(B[b][0]*B[b][0]+B[b][1]*B[b][1]+B[b][2]*B[b][2]));
      double n0=B[b][0]/L,n1=B[b][1]/L,n2=B[b][2]/L;
      for(int sg=-1;sg<=1;sg+=2){
        int ii=i+sg*B[b][0],jj=j+sg*B[b][1],kk=k+sg*B[b][2];
        if(ii<0||ii>=nx||jj<0||jj>=ny||kk<0||kk>=nz)continue;
        size_t qq=idx3(ii,jj,kk,ny,nz);
        double e=n0*(u[3*qq]-ux)+n1*(u[3*qq+1]-uy)+n2*(u[3*qq+2]-uz);
        f0+=e*n0; f1+=e*n1; f2+=e*n2;
      }
    }
    fu[3*pp]=f0;fu[3*pp+1]=f1;fu[3*pp+2]=f2;
  }
  const double r3=std::sqrt(3.0); int y0=ny/2,z0=nz/2;
  for(int i=0;i<nx-1;i++){
    double gi=g[i]; if(gi==0)continue;
    for(int oy=-1;oy<=0;oy++){double sy=oy==-1?-1.:1.;int j=y0+oy;
      for(int oz=-1;oz<=0;oz++){double sz=oz==-1?-1.:1.;int k=z0+oz;size_t p=idx3(i,j,k,ny,nz);
        double ux=u[3*p],uy=u[3*p+1],uz=u[3*p+2];
        double eL=(ux+sy*uy+sz*uz-q[i])/r3;
        double eR=(-ux+sy*uy+sz*uz+q[i+1])/r3;
        fu[3*p]+=-gi*(eL-eR)/r3; fu[3*p+1]+=-gi*sy*(eL+eR)/r3; fu[3*p+2]+=-gi*sz*(eL+eR)/r3;
        fq[i]+=gi*eL/r3; fq[i+1]+=-gi*eR/r3;
      }
    }
  }
}

double coupling_current(const std::vector<double>&q,const std::vector<double>&vq,const std::vector<double>&u,const std::vector<double>&vu,const std::vector<double>&g,int nx,int ny,int nz){
  const double r3=std::sqrt(3.0); int y0=ny/2,z0=nz/2; double J=0;
  for(int i=0;i<nx-1;i++){double gi=g[i];if(gi==0)continue;
    for(int oy=-1;oy<=0;oy++){double sy=oy==-1?-1.:1.;int j=y0+oy;
      for(int oz=-1;oz<=0;oz++){double sz=oz==-1?-1.:1.;int k=z0+oz;size_t p=idx3(i,j,k,ny,nz);
        double eL=(u[3*p]+sy*u[3*p+1]+sz*u[3*p+2]-q[i])/r3;
        double fbx=-gi*eL/r3,fby=-gi*eL*sy/r3,fbz=-gi*eL*sz/r3;
        J+=.5*(vu[3*p]*fbx+vu[3*p+1]*fby+vu[3*p+2]*fbz+vq[i]*fbx);
        double eR=(-u[3*p]+sy*u[3*p+1]+sz*u[3*p+2]+q[i+1])/r3;
        fbx=gi*eR/r3;fby=-gi*eR*sy/r3;fbz=-gi*eR*sz/r3;
        J+=.5*(vu[3*p]*fbx+vu[3*p+1]*fby+vu[3*p+2]*fbz+vq[i+1]*fbx);
      }
    }
  }
  return J;
}

Energies energies(const std::vector<double>&q,const std::vector<double>&vq,const std::vector<double>&u,const std::vector<double>&vu,const std::vector<double>&g,double Kq,int nx,int ny,int nz){
  Energies E; for(double x:vq)E.Eq+=.5*x*x; for(size_t p=0;p<vu.size();p++)E.Eb+=.5*vu[p]*vu[p];
  for(int i=0;i<nx-1;i++){double d=q[i+1]-q[i];E.Eq+=.5*Kq*d*d;}
  for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){
    size_t p=idx3(i,j,k,ny,nz);
    for(int b=0;b<NB;b++){int ii=i+B[b][0],jj=j+B[b][1],kk=k+B[b][2];if(ii<0||ii>=nx||jj<0||jj>=ny||kk<0||kk>=nz)continue;
      size_t qq=idx3(ii,jj,kk,ny,nz);double L=std::sqrt(double(B[b][0]*B[b][0]+B[b][1]*B[b][1]+B[b][2]*B[b][2]));
      double e=(B[b][0]*(u[3*qq]-u[3*p])+B[b][1]*(u[3*qq+1]-u[3*p+1])+B[b][2]*(u[3*qq+2]-u[3*p+2]))/L;E.Eb+=.5*e*e;
    }
  }
  const double r3=std::sqrt(3.0);int y0=ny/2,z0=nz/2;
  for(int i=0;i<nx-1;i++){double gi=g[i];if(gi==0)continue;for(int oy=-1;oy<=0;oy++){double sy=oy==-1?-1.:1.;int j=y0+oy;for(int oz=-1;oz<=0;oz++){double sz=oz==-1?-1.:1.;int k=z0+oz;size_t p=idx3(i,j,k,ny,nz);
    double eL=(u[3*p]+sy*u[3*p+1]+sz*u[3*p+2]-q[i])/r3;double eR=(-u[3*p]+sy*u[3*p+1]+sz*u[3*p+2]+q[i+1])/r3;E.Ec+=.5*gi*(eL*eL+eR*eR);
  }}}
  return E;
}

double historical_shell(const std::vector<double>&vu,int nx,int ny,int nz,int w){double den=0,num=0;for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){size_t p=idx3(i,j,k,ny,nz);double e=vu[3*p]*vu[3*p]+vu[3*p+1]*vu[3*p+1]+vu[3*p+2]*vu[3*p+2];den+=e;if(j<w||j>=ny-w||k<w||k>=nz-w)num+=e;}return den?num/den:0;}
Shell shell_total(const std::vector<double>&u,const std::vector<double>&vu,int nx,int ny,int nz,int w){
  size_t nc=(size_t)nx*ny*nz;std::vector<double> en(nc,0.0);
  for(size_t p=0;p<nc;p++)en[p]=.5*(vu[3*p]*vu[3*p]+vu[3*p+1]*vu[3*p+1]+vu[3*p+2]*vu[3*p+2]);
  for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){size_t p=idx3(i,j,k,ny,nz);for(int b=0;b<NB;b++){int ii=i+B[b][0],jj=j+B[b][1],kk=k+B[b][2];if(ii<0||ii>=nx||jj<0||jj>=ny||kk<0||kk>=nz)continue;size_t qq=idx3(ii,jj,kk,ny,nz);double L=std::sqrt(double(B[b][0]*B[b][0]+B[b][1]*B[b][1]+B[b][2]*B[b][2]));double e=(B[b][0]*(u[3*qq]-u[3*p])+B[b][1]*(u[3*qq+1]-u[3*p+1])+B[b][2]*(u[3*qq+2]-u[3*p+2]))/L;double h=.25*e*e;en[p]+=h;en[qq]+=h;}}
  Shell s;for(int i=0;i<nx;i++)for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){size_t p=idx3(i,j,k,ny,nz);double e=en[p];s.node_sum+=e;std::array<bool,6> f={i<w,i>=nx-w,j<w,j>=ny-w,k<w,k>=nz-w};bool any=false;for(int z=0;z<6;z++){if(f[z]){s.face[z]+=e;any=true;}}if(any)s.union_total+=e;}return s;
}

struct CheckHeader{
  uint64_t magic=0x50303743484b5031ULL; // P07CHKP1
  int32_t nx=0,ny=0,nz=0;
  int64_t step=0;
  double dt=0,Etot0=0,E0=0,cum=0,Jprev=0;
};

template<class T> void write_vec(std::ofstream&f,const std::vector<T>&v){uint64_t n=v.size();f.write((char*)&n,8);f.write((char*)v.data(),n*sizeof(T));}
template<class T> void read_vec(std::ifstream&f,std::vector<T>&v){uint64_t n=0;f.read((char*)&n,8);v.resize(n);f.read((char*)v.data(),n*sizeof(T));}

void save_checkpoint(const std::string&path,const CheckHeader&h,const std::vector<double>&q,const std::vector<double>&vq,const std::vector<double>&u,const std::vector<double>&vu){
  std::ofstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot open checkpoint output");f.write((char*)&h,sizeof(h));write_vec(f,q);write_vec(f,vq);write_vec(f,u);write_vec(f,vu);if(!f)throw std::runtime_error("checkpoint write failed");
}
void load_checkpoint(const std::string&path,CheckHeader&h,std::vector<double>&q,std::vector<double>&vq,std::vector<double>&u,std::vector<double>&vu){
  std::ifstream f(path,std::ios::binary);if(!f)throw std::runtime_error("cannot open checkpoint input");f.read((char*)&h,sizeof(h));if(h.magic!=0x50303743484b5031ULL)throw std::runtime_error("bad checkpoint magic");read_vec(f,q);read_vec(f,vq);read_vec(f,u);read_vec(f,vu);if(!f)throw std::runtime_error("checkpoint read failed");
}

int main(int argc,char**argv){
  if(argc<13){std::cerr<<"usage: "<<argv[0]<<" THREADS CHECKPOINT_IN|- CHECKPOINT_OUT|- OUTPUT_CSV|- OUTPUT_BIN|- NX NY NZ DT END_T SAMPLE_START SAMPLE_STEP\n";return 2;}
  int threads=std::atoi(argv[1]);omp_set_num_threads(threads);std::string chk_in=argv[2],chk_out=argv[3],csv_path=argv[4],bin_path=argv[5];
  int nx=std::atoi(argv[6]),ny=std::atoi(argv[7]),nz=std::atoi(argv[8]);double dt=std::atof(argv[9]),end_t=std::atof(argv[10]),sample_start=std::atof(argv[11]),sample_step=std::atof(argv[12]);
  const double Kq=1.69;const int end_step=std::llround(end_t/dt);const int stride=std::llround(sample_step/dt);if(stride<=0||std::abs(stride*dt-sample_step)>1e-12)return 3;
  std::vector<double>q,vq,u,vu;CheckHeader h;size_t nv=(size_t)nx*ny*nz*3;std::vector<double>g(nx-1,0.0),fq(nx),fu(nv);double gs[6]={.1,.2,.2,.2,.2,.1};for(int s=0;s<6;s++)g[60+s]=gs[s];
  if(chk_in!="-"){
    load_checkpoint(chk_in,h,q,vq,u,vu);if(h.nx!=nx||h.ny!=ny||h.nz!=nz||std::abs(h.dt-dt)>0)throw std::runtime_error("checkpoint geometry/dt mismatch");
  }else{
    h.nx=nx;h.ny=ny;h.nz=nz;h.dt=dt;h.step=0;q.assign(nx,0);vq.assign(nx,0);u.assign(nv,0);vu.assign(nv,0);
    const double k0=2*std::asin(1.2/(2*std::sqrt(Kq)));for(int i=0;i<nx;i++)q[i]=.002*std::exp(-.5*std::pow((i-18.)/7.,2))*std::cos(k0*(i-18.));
    const std::complex<double>I(0,1);std::vector<std::complex<double>>qh(nx),vh(nx);for(int m=0;m<nx;m++){double kk=2*M_PI*(m<=nx/2?m:m-nx)/nx;for(int n=0;n<nx;n++)qh[m]+=q[n]*std::exp(-I*kk*double(n));double om=2*std::sqrt(Kq)*std::abs(std::sin(kk/2));double sg=(m==0||m==nx/2)?0:(kk>0?1:-1);vh[m]=-I*sg*om*qh[m];}for(int n=0;n<nx;n++){std::complex<double>z=0;for(int m=0;m<nx;m++){double kk=2*M_PI*(m<=nx/2?m:m-nx)/nx;z+=vh[m]*std::exp(I*kk*double(n));}vq[n]=z.real()/nx;}
    force_all(q,u,g,Kq,nx,ny,nz,fq,fu);Energies Einit=energies(q,vq,u,vu,g,Kq,nx,ny,nz);h.Etot0=Einit.Eq+Einit.Eb+Einit.Ec;h.E0=Einit.Eq;h.Jprev=coupling_current(q,vq,u,vu,g,nx,ny,nz);h.cum=0;
  }
  force_all(q,u,g,Kq,nx,ny,nz,fq,fu);
  bool write_output=(csv_path!="-"&&bin_path!="-");std::ofstream csv;std::ofstream bin;std::vector<float>tmp;int written=0;
  if(write_output){csv.open(csv_path);csv<<"time,T_current,T_coupling,T_bulk_open,current,drift,historical_shell,total_shell_union,total_shell_xlow,total_shell_xhigh,total_shell_ylow,total_shell_yhigh,total_shell_zlow,total_shell_zhigh,node_closure\n"<<std::setprecision(17);bin.open(bin_path,std::ios::binary);uint64_t ns=(uint64_t)std::llround((end_t-sample_start)/sample_step)+1,a=nx,b=ny,c=nz;bin.write((char*)&ns,8);bin.write((char*)&a,8);bin.write((char*)&b,8);bin.write((char*)&c,8);tmp.resize(nv);}
  double start=omp_get_wtime();
  for(int64_t st=h.step+1;st<=end_step;st++){
    #pragma omp parallel for schedule(static)
    for(int i=0;i<nx;i++)vq[i]+=.5*dt*fq[i];
    #pragma omp parallel for schedule(static)
    for(long long i=0;i<(long long)nv;i++)vu[i]+=.5*dt*fu[i];
    #pragma omp parallel for schedule(static)
    for(int i=0;i<nx;i++)q[i]+=dt*vq[i];
    #pragma omp parallel for schedule(static)
    for(long long i=0;i<(long long)nv;i++)u[i]+=dt*vu[i];
    force_all(q,u,g,Kq,nx,ny,nz,fq,fu);
    #pragma omp parallel for schedule(static)
    for(int i=0;i<nx;i++)vq[i]+=.5*dt*fq[i];
    #pragma omp parallel for schedule(static)
    for(long long i=0;i<(long long)nv;i++)vu[i]+=.5*dt*fu[i];
    double J=coupling_current(q,vq,u,vu,g,nx,ny,nz);h.cum+=.5*dt*(h.Jprev+J);h.Jprev=J;h.step=st;
    double t=st*dt;if(write_output&&t>=sample_start-1e-12&&st%stride==0){Energies E=energies(q,vq,u,vu,g,Kq,nx,ny,nz);double hs=historical_shell(vu,nx,ny,nz,2);Shell sh;bool integer=std::abs(t-std::round(t))<1e-10&&t>=52-1e-10&&t<=55+1e-10;if(integer)sh=shell_total(u,vu,nx,ny,nz,2);double nodeclose=integer?(sh.node_sum-E.Eb)/E.Eb:0;csv<<t<<","<<h.cum/h.E0<<","<<E.Ec/h.E0<<","<<E.Eb/h.E0<<","<<J/h.E0<<","<<(E.Eq+E.Eb+E.Ec-h.Etot0)/h.Etot0<<","<<hs<<","<<(integer?sh.union_total/E.Eb:-1);for(int z=0;z<6;z++)csv<<","<<(integer?sh.face[z]/E.Eb:-1);csv<<","<<nodeclose<<"\n";bin.write((char*)&t,sizeof(double));for(size_t i=0;i<nv;i++)tmp[i]=(float)u[i];bin.write((char*)tmp.data(),tmp.size()*sizeof(float));for(size_t i=0;i<nv;i++)tmp[i]=(float)vu[i];bin.write((char*)tmp.data(),tmp.size()*sizeof(float));written++;std::cerr<<"sample "<<t<<" Eb "<<E.Eb/h.E0<<" shell "<<(integer?sh.union_total/E.Eb:-1)<<"\n";}
  }
  if(chk_out!="-")save_checkpoint(chk_out,h,q,vq,u,vu);std::cerr<<"segment runtime "<<(omp_get_wtime()-start)<<" end_step "<<h.step<<" written "<<written<<"\n";return 0;
}
