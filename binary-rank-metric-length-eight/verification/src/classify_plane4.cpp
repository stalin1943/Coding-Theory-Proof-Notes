// Exhaustively classify all rank-four binary subspaces of F64^2 with
// point weights at most two, up to GL(2,F64). Every normalized member
// contains (1,0),(0,1); equivalent normalized members are explicitly marked.
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
int M[64][64],inv_[64],pid[4096];
int mul(int a,int b){int c=0;while(b){if(b&1)c^=a;b>>=1;a<<=1;if(a&64)a^=91;}return c;}
int lift(int x){return ((x&31)<<1)|((x>>5)<<7);}
int project(int x){return ((x>>1)&31)|((x>>7)<<5);}
int main(){
 for(int a=0;a<64;++a)for(int b=0;b<64;++b)M[a][b]=mul(a,b);
 for(int a=1;a<64;++a)for(int b=1;b<64;++b)if(M[a][b]==1)inv_[a]=b;
 for(int v=1;v<4096;++v){int x=v&63,y=v>>6;pid[v]=x?M[y][inv_[x]]:64;}
 vector<unsigned char>seen(1<<20),valid(1<<20);long long total=0,good=0;map<int,int>hist;vector<pair<int,int>>all;
 for(int a=1;a<1024;++a)for(int b=a+1;b<1024;++b)if(b<(a^b)){
  ++total;int cnt[65]={};int bas[4]={1,64,lift(a),lift(b)};vector<int>T{0};for(int v:bas){int n=T.size();for(int i=0;i<n;++i)T.push_back(T[i]^v);}bool ok=true;int s=0;for(int v:T)if(v){int n=++cnt[pid[v]];if(n==4){ok=false;break;}if(n==2)++s;}
  if(ok){++good;valid[(a<<10)|b]=s+1;all.push_back({a,b});++hist[s];}
 }
 ofstream cat("plane4_catalog.csv"),cover("plane4_cover.csv");cat<<"orbit,a,b,point_weight2,normalized_members,stabilizer_order\n";cover<<"orbit,a,b\n";map<int,int>orbhist;int orb=0,marked=0;
 for(auto [a,b]:all){int key=(a<<10)|b;if(seen[key])continue;vector<int>T{0};int bas[4]={1,64,lift(a),lift(b)};for(int v:bas){int n=T.size();for(int i=0;i<n;++i)T.push_back(T[i]^v);}vector<int>keys;int frames=0;
  for(int u:T)if(u)for(int v:T)if(v){int det=M[u&63][v>>6]^M[u>>6][v&63];if(!det)continue;++frames;int id=inv_[det],x0=0,y0=0;
   for(int w:bas){int x=M[id][M[v>>6][w&63]^M[v&63][w>>6]];int y=M[id][M[u>>6][w&63]^M[u&63][w>>6]];int q=project(x|(y<<6));if(!q)continue;if(!x0)x0=q;else if(q!=x0&&!y0)y0=q;}
   assert(x0&&y0);array<int,3>q{x0,y0,x0^y0};sort(q.begin(),q.end());int kk=(q[0]<<10)|q[1];assert(valid[kk]==valid[key]);keys.push_back(kk);
  }
  sort(keys.begin(),keys.end());keys.erase(unique(keys.begin(),keys.end()),keys.end());assert(!keys.empty()&&keys.front()==key);assert(frames%keys.size()==0);
  cat<<orb<<','<<a<<','<<b<<','<<(int)valid[key]-1<<','<<keys.size()<<','<<frames/keys.size()<<'\n';++orbhist[valid[key]-1];
  for(int kk:keys){assert(!seen[kk]);seen[kk]=1;++marked;cover<<orb<<','<<(kk>>10)<<','<<(kk&1023)<<'\n';}++orb;
 }
 assert(marked==good&&total==174251);
 cout<<"normalized_total "<<total<<" good "<<good<<" orbits "<<orb<<" covered "<<marked<<'\n';for(auto [s,n]:hist)cout<<"heavy "<<s<<" normalized "<<n<<" orbits "<<orbhist[s]<<'\n';
}
