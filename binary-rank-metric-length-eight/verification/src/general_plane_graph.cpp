// Complete 8D graph lifts over normalized 4D hyperplane intersections.
// Assumptions: U subset F64^3 has binary dimension 8, all F-hyperplane
// intersections have dimension <=4, all point intersections have dim <=2.
// The independently generated plane4_catalog covers every possible T=U cap H.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
int M[64][64],inv_[64],pid[1<<18],iso[64],poly=91;
using u64=uint64_t;
struct Bin{
 int n;vector<u64>row;vector<int>free;
 Bin(int n_=0):n(n_),row(n_,0){}
 bool add(u64 x){for(int j=n-1;j>=0;--j)if((x>>j)&1){if(row[j])x^=row[j];else{row[j]=x;return true;}}return false;}
 u64 rem(u64 x)const{for(int j=n-1;j>=0;--j)if(((x>>j)&1)&&row[j])x^=row[j];return x;}
 void finish(){free.clear();for(int j=0;j<n;++j)if(!row[j])free.push_back(j);}
 u64 pack(u64 x)const{x=rem(x);u64 y=0;for(int k=0;k<(int)free.size();++k)y|=((x>>free[k])&1)<<k;return y;}
 u64 unpack(u64 x)const{u64 y=0;for(int k=0;k<(int)free.size();++k)y|=((x>>k)&1)<<free[k];return y;}
 int rank()const{int r=0;for(u64 x:row)r+=(x!=0);return r;}
};
vector<int> span(const vector<int>&b){vector<int>s{0};for(int x:b){int n=s.size();for(int i=0;i<n;++i)s.push_back(s[i]^x);}return s;}
int mulraw(int a,int b,int m){int r=0;while(b){if(b&1)r^=a;b>>=1;a<<=1;if(a&64)a^=m;}return r;}
int dot(int h,int v){return M[h&63][v&63]^M[(h>>6)&63][(v>>6)&63]^M[h>>12][v>>12];}
int cross_(int a,int b){return (M[(a>>6)&63][b>>12]^M[a>>12][(b>>6)&63])|((M[a>>12][b&63]^M[a&63][b>>12])<<6)|((M[a&63][(b>>6)&63]^M[(a>>6)&63][b&63])<<12);}
int trace_(int a){int t=0;for(int i=0;i<6;++i){t^=a;a=M[a][a];}assert(t==0||t==1);return t;}
int lift(int x){return ((x&31)<<1)|((x>>5)<<7);}
int map12(int x){return iso[x&63]|(iso[x>>6]<<6);}
using Map=array<int,4>;
pair<vector<Map>,vector<Map>> graphmaps(const vector<int>&T,const vector<int>&Z){
 Bin qt(12);for(int x:T)qt.add(x);qt.finish();assert(qt.free.size()==8);
 Bin sh(32);for(int c=0;c<2;++c)for(int k=0;k<6;++k){u64 code=0;for(int j=0;j<4;++j)code|=qt.pack(M[1<<k][Z[j]]<<(6*c))<<(8*j);sh.add(code);}sh.finish();assert(sh.rank()==12&&sh.free.size()==20);
 vector<Map>gen(20);for(int k=0;k<20;++k){u64 code=sh.unpack(u64(1)<<k);for(int j=0;j<4;++j)gen[k][j]=qt.unpack((code>>(8*j))&255);}
 vector<Map>left(1024),right(1024);for(int x=1;x<1024;++x){int k=__builtin_ctz((unsigned)x),y=x&(x-1);for(int j=0;j<4;++j){left[x][j]=left[y][j]^gen[k][j];right[x][j]=right[y][j]^gen[k+10][j];}}
 return {left,right};
}
struct Dir{int a,b,r,d;vector<int>sx,sy;vector<unsigned char>bad;};
Dir direction(int a,int b,const vector<int>&T,const vector<int>&Z,const vector<Map>&left,const vector<Map>&right){
 auto ell=[&](int v){return M[a][v&63]^M[b][v>>6];};
 Bin S(6);for(int v:T)S.add(ell(v));S.finish();int c=S.free.size();assert(c>=2&&c<=4);
 Bin R(4*c);for(int k=0;k<6;++k){int code=0;for(int j=0;j<4;++j)code|=S.pack(M[1<<k][Z[j]])<<(c*j);R.add(code);}R.finish();
 Dir d{a,b,S.rank(),(int)R.free.size(),vector<int>(1024),vector<int>(1024),{}};
 for(int x=0;x<1024;++x){int lc=0,rc=0;for(int j=0;j<4;++j){lc|=S.pack(ell(left[x][j]))<<(c*j);rc|=S.pack(ell(right[x][j]))<<(c*j);}d.sx[x]=R.pack(lc);d.sy[x]=R.pack(rc);}
 d.bad.assign(1<<d.d,0);d.bad[0]=1;if(d.r==2){assert(c==4);for(int u=1;u<16;++u)for(int v=1;v<16;++v){int mat=0;for(int j=0;j<4;++j)if((v>>j)&1)mat|=u<<(4*j);d.bad[R.pack(mat)]=1;}}
 return d;
}
vector<int>points;
// Direct finite-field rank and incidence verification. No syndrome tables.
int direct(const vector<int>&b,int&heavy,int&h0,int&h1){
 auto U=span(b);vector<int>count(4161);for(int x:U)if(x)++count[pid[x]];heavy=0;for(int n:count){heavy+=(n==3);if(n>3)return 1;}
 vector<int>h4;for(int h:points){Bin im(6);for(int x:b)im.add(dot(h,x));if(im.rank()<4){h0=h;h1=0;return 2;}if(im.rank()==4)h4.push_back(h);}
 assert((int)h4.size()==31+4*heavy);
 for(int i=0;i<(int)h4.size();++i)for(int j=i+1;j<(int)h4.size();++j)if(!count[pid[cross_(h4[i],h4[j])]]){h0=h4[i];h1=h4[j];return 3;}
 return 0;
}
int main(int argc,char**argv){
 int lo=argc>1?stoi(argv[1]):0,hi=argc>2?stoi(argv[2]):911;poly=argc>3?stoi(argv[3]):91;
 for(int a=0;a<64;++a)for(int b=0;b<64;++b)M[a][b]=mulraw(a,b,poly);for(int a=1;a<64;++a)for(int b=1;b<64;++b)if(M[a][b]==1)inv_[a]=b;
 for(int a=0;a<64;++a){int z=0;for(int j=5;j>=0;--j)z=M[z][poly==91?2:7]^((a>>j)&1);iso[a]=z;}
 for(int a=0;a<64;++a)for(int b=0;b<64;++b)assert(iso[mulraw(a,b,91)]==M[iso[a]][iso[b]]);
 int bs0[11]={14,2,4,6,8,12,16,20,22,28,36},bs[11];for(int i=0;i<11;++i)bs[i]=iso[bs0[i]];
 for(int y=0;y<64;++y)for(int z=0;z<64;++z)points.push_back(1|(y<<6)|(z<<12));for(int z=0;z<64;++z)points.push_back(64|(z<<12));points.push_back(4096);
 for(int v=1;v<(1<<18);++v){int x=v&63,y=(v>>6)&63,z=v>>12;pid[v]=x?64*M[y][inv_[x]]+M[z][inv_[x]]:y?4096+M[z][inv_[y]]:4160;}
 vector<vector<int>>Zs;for(int k=0;k<11;++k){vector<int>Z;Bin s(6);for(int x=1;x<64;++x)if(!trace_(x)&&!trace_(M[bs[k]][x]))if(s.add(x))Z.push_back(x);assert(Z.size()==4);Zs.push_back(Z);}
 string tag=to_string(lo)+"_"+to_string(hi)+"_"+to_string(poly);ofstream out("general_counts_"+tag+".csv"),cand("general_survivors_"+tag+".txt");
 out<<"orbit,k,core_heavy,normal_forms,after_rank34,after_rank2,point_bad,hyperplane_bad,two_span,target\n";
 ifstream fcat("plane4_catalog.csv");string line;getline(fcat,line);long long all=0,n1all=0,n2all=0,targets=0;int done=0;auto start=chrono::steady_clock::now();
 while(getline(fcat,line)){replace(line.begin(),line.end(),',',' ');istringstream sin(line);int orb,a,b,th;sin>>orb>>a>>b>>th;if(orb<lo||orb>hi)continue;
  vector<int>T{1,64,map12(lift(a)),map12(lift(b))};
  for(int k=0;k<11;++k){auto&Z=Zs[k];auto [left,right]=graphmaps(T,Z);vector<Dir>ds;
   for(int d=0;d<64;++d)ds.push_back(direction(1,d,T,Z,left,right));ds.push_back(direction(0,1,T,Z,left,right));stable_sort(ds.begin(),ds.end(),[](const Dir&a,const Dir&b){return a.r>b.r;});
   array<array<u64,16>,1024>allowed;for(auto&row:allowed)row.fill(~u64(0));vector<int>active(1024);iota(active.begin(),active.end(),0);
   for(auto&d:ds){if(d.r==2)continue;vector<array<u64,16>>group(1<<d.d);for(int y=0;y<1024;++y)group[d.sy[y]][y>>6]|=u64(1)<<(y&63);int n=0;
    for(int x:active){bool nz=false;auto&g=group[d.sx[x]];for(int w=0;w<16;++w){allowed[x][w]&=~g[w];nz|=(allowed[x][w]!=0);}if(nz)active[n++]=x;}active.resize(n);if(active.empty())break;
   }
   long long n1=0,n2=0,bad[4]={};for(int x:active)for(int w=0;w<16;++w){u64 bits=allowed[x][w];while(bits){int y=w*64+__builtin_ctzll(bits);bits&=bits-1;++n1;bool fail=false;
     for(auto&d:ds)if(d.r==2&&d.bad[d.sx[x]^d.sy[y]]){fail=true;break;}if(fail)continue;++n2;
     vector<int>basis=T;for(int j=0;j<4;++j)basis.push_back((left[x][j]^right[y][j])|(Z[j]<<12));int h=0,h0=0,h1=0;int result=direct(basis,h,h0,h1);++bad[result];
     cand<<orb<<' '<<k<<' '<<th<<' '<<x<<' '<<y<<' '<<result<<' '<<h<<' '<<h0<<' '<<h1;for(int v:basis)cand<<' '<<v;cand<<'\n';
   }}
   out<<orb<<','<<k<<','<<th<<",1048576,"<<n1<<','<<n2<<','<<bad[1]<<','<<bad[2]<<','<<bad[3]<<','<<bad[0]<<'\n';out.flush();cand.flush();all+=1048576;n1all+=n1;n2all+=n2;targets+=bad[0];++done;
   if(bad[0])cerr<<"TARGET at orbit "<<orb<<" k "<<k<<" n "<<bad[0]<<'\n';
  }
  if(orb%20==0)cerr<<"orbit "<<orb<<" cases "<<done<<" stage1 "<<n1all<<" stage2 "<<n2all<<" targets "<<targets<<" elapsed "<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<'\n';
 }
 cout<<"cases "<<done<<" normal_forms "<<all<<" stage1 "<<n1all<<" stage2 "<<n2all<<" targets "<<targets<<'\n';
}
