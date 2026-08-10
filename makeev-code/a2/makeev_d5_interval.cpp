/*
  Certified real branch-and-bound test for the proposed d=5 Makeev
  counterexample.  The input polynomials encode the multiplication tensor of
  a normalized six-point simplex frame.  A successful exhaustive run prints

      Certificate verified: no real tensor zero exists.

  and exits 0.  Exit 2 means only that the requested box/depth budget was
  exhausted.  Compile without fast-math; every elementary floating-point
  operation is enlarged by one ulp in each direction.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

#include "tensor_polynomials.inc"

static_assert(std::numeric_limits<double>::is_iec559,
              "This verifier requires IEEE-754 binary floating point.");

namespace {

constexpr int RK=20;

double down(double x) {
#ifdef FAST_PROTO
  return x;
#else
  return std::nextafter(x,-std::numeric_limits<double>::infinity());
#endif
}
double up(double x) {
#ifdef FAST_PROTO
  return x;
#else
  return std::nextafter(x,std::numeric_limits<double>::infinity());
#endif
}
struct I {
  double lo=0,hi=0;
  I()=default;
  explicit I(double x):lo(down(x)),hi(up(x)){}
  I(double a,double b):lo(a),hi(b){}
  bool zero()const{return lo<=0&&hi>=0;}
  bool empty()const{return lo>hi;}
  double mid()const{return .5*lo+.5*hi;}
  double width()const{return hi-lo;}
};
I operator+(I a,I b){return {down(a.lo+b.lo),up(a.hi+b.hi)};}
I operator-(I a,I b){return {down(a.lo-b.hi),up(a.hi-b.lo)};}
I operator*(I a,I b){double p[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};return{down(*std::min_element(p,p+4)),up(*std::max_element(p,p+4))};}
I recip(I a){if(a.zero())return{-INFINITY,INFINITY};return{down(1/a.hi),up(1/a.lo)};}
I operator/(I a,I b){return a*recip(b);}
I cap(I a,I b){return{std::max(a.lo,b.lo),std::min(a.hi,b.hi)};}

struct Eval {
  I value;
  I pcenter;
  std::array<I,TP_NV> grad{};
  double pmid=0;
  std::array<double,TP_NV> gmid{};
};

int mono_key(const unsigned char* v,int d){
  int key=d;
  for(int i=0;i<d;i++)key|=(int(v[i])+1)<<(2+5*i);
  return key;
}
I exact_taylor_combination(const std::array<double,TP_NE>& weight,
                           const std::array<double,TP_NV>& mid,
                           const std::array<I,TP_NV>& dx) {
  std::unordered_map<int,I> shifted;
  for(int e=0;e<TP_NE;e++)if(weight[e]!=0.0)
    for(int ti=TP_OFF[e];ti<TP_OFF[e+1];ti++){
      const Term&t=TP_TERMS[ti];unsigned char v[3]={t.v0,t.v1,t.v2};
      for(int mask=0;mask<(1<<t.degree);mask++){
        I c=I(weight[e])*I(double(t.c));unsigned char chosen[3]{};int d=0;
        for(int k=0;k<t.degree;k++)if(mask&(1<<k))chosen[d++]=v[k];else c=c*I(mid[v[k]]);
        int key=mono_key(chosen,d);shifted[key]=shifted[key]+c;
      }
    }
  I ans;
  for(auto [key,c]:shifted){int d=key&3;I term=c;
    for(int k=0;k<d;k++){int v=((key>>(2+5*k))&31)-1;term=term*dx[v];}
    ans=ans+term;
  }
  return ans;
}

I exact_taylor_poly(int e,const std::array<double,TP_NV>& mid,
                    const std::array<I,TP_NV>& dx) {
  std::unordered_map<int,I> shifted;
  for(int ti=TP_OFF[e];ti<TP_OFF[e+1];ti++){
    const Term&t=TP_TERMS[ti];unsigned char v[3]={t.v0,t.v1,t.v2};
    for(int mask=0;mask<(1<<t.degree);mask++){
      I c(double(t.c));unsigned char chosen[3]{};int d=0;
      for(int k=0;k<t.degree;k++)if(mask&(1<<k))chosen[d++]=v[k];else c=c*I(mid[v[k]]);
      int key=mono_key(chosen,d);shifted[key]=shifted[key]+c;
    }
  }
  I ans;
  for(auto [key,c]:shifted){int d=key&3;I term=c;
    for(int k=0;k<d;k++){int v=((key>>(2+5*k))&31)-1;term=term*dx[v];}
    ans=ans+term;
  }
  return ans;
}

Eval eval_poly(int poly,const std::array<I,TP_NV>& x,
               const std::array<double,TP_NV>& m) {
  Eval e;
  for(int ti=TP_OFF[poly];ti<TP_OFF[poly+1];++ti) {
    const Term &t=TP_TERMS[ti];
    I iv(double(t.c)),pc(double(t.c)); double pv=double(t.c);
    unsigned char vv[3]={t.v0,t.v1,t.v2};
    for(int k=0;k<t.degree;k++){iv=iv*x[vv[k]];pc=pc*I(m[vv[k]]);pv*=m[vv[k]];}
    e.value=e.value+iv;e.pcenter=e.pcenter+pc;e.pmid+=pv;
    for(int k=0;k<t.degree;k++) {
      I di(double(t.c));double dp=double(t.c);
      for(int l=0;l<t.degree;l++)if(l!=k){di=di*x[vv[l]];dp*=m[vv[l]];}
      e.grad[vv[k]]=e.grad[vv[k]]+di;e.gmid[vv[k]]+=dp;
    }
  }
  return e;
}

bool inverse(double a[RK][RK],double z[RK][RK]) {
  double q[RK][2*RK]{};
  for(int i=0;i<RK;i++)for(int j=0;j<RK;j++){q[i][j]=a[i][j];q[i][j+RK]=(i==j);}
  for(int k=0;k<RK;k++){
    int p=k;for(int i=k+1;i<RK;i++)if(std::abs(q[i][k])>std::abs(q[p][k]))p=i;
    if(std::abs(q[p][k])<1e-11)return false;
    if(p!=k)for(int j=0;j<2*RK;j++)std::swap(q[p][j],q[k][j]);
    double d=q[k][k];for(int j=0;j<2*RK;j++)q[k][j]/=d;
    for(int i=0;i<RK;i++)if(i!=k){double s=q[i][k];for(int j=0;j<2*RK;j++)q[i][j]-=s*q[k][j];}
  }
  for(int i=0;i<RK;i++)for(int j=0;j<RK;j++)z[i][j]=q[i][j+RK];
  return true;
}

std::array<int,RK> select_rows(const std::array<Eval,TP_NE>& f) {
  std::array<int,RK> ans{};std::array<bool,TP_NE> used{};
  std::array<std::array<double,TP_NV>,RK> basis{};
  for(int k=0;k<RK;k++){
    int best=-1;double bn=-1;
    for(int i=0;i<TP_NE;i++)if(!used[i]){
      auto w=f[i].gmid;
      for(int h=0;h<k;h++){double d=0;for(int j=0;j<TP_NV;j++)d+=w[j]*basis[h][j];for(int j=0;j<TP_NV;j++)w[j]-=d*basis[h][j];}
      double n=0;for(double y:w)n+=y*y;if(n>bn){bn=n;best=i;}
    }
    ans[k]=best;used[best]=true;auto w=f[best].gmid;
    for(int h=0;h<k;h++){double d=0;for(int j=0;j<TP_NV;j++)d+=w[j]*basis[h][j];for(int j=0;j<TP_NV;j++)w[j]-=d*basis[h][j];}
    double n=0;for(double y:w)n+=y*y;n=std::sqrt(n);if(n)for(int j=0;j<TP_NV;j++)basis[k][j]=w[j]/n;
  }
  return ans;
}

std::array<int,RK> select_cols(const std::array<Eval,TP_NE>& f,
                               const std::array<int,RK>& rows) {
  std::array<int,RK> ans{};std::array<bool,TP_NV> used{};
  std::array<std::array<double,RK>,RK> basis{};
  for(int k=0;k<RK;k++){
    int best=-1;double bn=-1;
    for(int col=0;col<TP_NV;col++)if(!used[col]){
      std::array<double,RK>w{};for(int i=0;i<RK;i++)w[i]=f[rows[i]].gmid[col];
      for(int h=0;h<k;h++){double d=0;for(int i=0;i<RK;i++)d+=w[i]*basis[h][i];for(int i=0;i<RK;i++)w[i]-=d*basis[h][i];}
      double n=0;for(double y:w)n+=y*y;if(n>bn){bn=n;best=col;}
    }
    ans[k]=best;used[best]=true;std::array<double,RK>w{};
    for(int i=0;i<RK;i++)w[i]=f[rows[i]].gmid[best];
    for(int h=0;h<k;h++){double d=0;for(int i=0;i<RK;i++)d+=w[i]*basis[h][i];for(int i=0;i<RK;i++)w[i]-=d*basis[h][i];}
    double n=0;for(double y:w)n+=y*y;n=std::sqrt(n);if(n)for(int i=0;i<RK;i++)basis[k][i]=w[i]/n;
  }
  return ans;
}

struct Box{std::array<I,TP_NV>x;int depth=0;};
enum R{DROP,KEEP};
R process(Box& b,std::uint64_t& range,std::uint64_t& lin,
          std::uint64_t& invfail,std::uint64_t& nw,int& split_hint) {
  const double tb=1.0;const I allowed(-tb,tb);
  // Every cubic moment is bounded by sqrt(5/6).  Apply the bound not only
  // to the 21 free moments but also to the 14 moments eliminated linearly.
  for(int pass=0;pass<2;pass++)for(int e=0;e<35;e++){
    I total;
    for(int k=TP_LOFF[e];k<TP_LOFF[e+1];k++){
      const auto&t=TP_LIN[k];total=total+I(double(t.num)/t.den)*b.x[t.var];
    }
    if(cap(total,allowed).empty()){range++;return DROP;}
    for(int k=TP_LOFF[e];k<TP_LOFF[e+1];k++){
      const auto&t=TP_LIN[k];I rest;
      for(int l=TP_LOFF[e];l<TP_LOFF[e+1];l++)if(l!=k){const auto&s=TP_LIN[l];rest=rest+I(double(s.num)/s.den)*b.x[s.var];}
      I take=(allowed-rest)/I(double(t.num)/t.den);
      b.x[t.var]=cap(b.x[t.var],take);if(b.x[t.var].empty()){range++;return DROP;}
    }
  }
  std::array<double,TP_NV> m{};std::array<I,TP_NV> dx{};
  for(int i=0;i<TP_NV;i++){m[i]=b.x[i].mid();dx[i]=b.x[i]-I(m[i]);}
  std::array<Eval,TP_NE> f{};
  for(int i=0;i<TP_NE;i++){
    f[i]=eval_poly(i,b.x,m);I centered=exact_taylor_poly(i,m,dx);
    f[i].value=cap(f[i].value,centered);
    if(f[i].value.empty()||!f[i].value.zero()){range++;return DROP;}
  }
  // Stable overdetermined consistency test.  Project F(mid) orthogonally
  // away from the column space of J(mid); the normalized residual is a
  // bounded left-null multiplier, unlike multipliers obtained from a
  // potentially ill-conditioned square minor.
  std::array<std::array<double,TP_NE>,RK> qbasis{};
  std::array<bool,TP_NV> used_col{};int qrank=0;
  for(int step=0;step<RK;step++){
    int best=-1;double bn=-1;std::array<double,TP_NE> bw{};
    for(int col=0;col<TP_NV;col++)if(!used_col[col]){
      std::array<double,TP_NE>w{};for(int i=0;i<TP_NE;i++)w[i]=f[i].gmid[col];
      for(int h=0;h<qrank;h++){double d=0;for(int i=0;i<TP_NE;i++)d+=w[i]*qbasis[h][i];for(int i=0;i<TP_NE;i++)w[i]-=d*qbasis[h][i];}
      double n=0;for(double y:w)n+=y*y;if(n>bn){bn=n;best=col;bw=w;}
    }
    if(best<0||bn<1e-22)break;
    used_col[best]=true;
    double n=std::sqrt(bn);
    for(int i=0;i<TP_NE;i++)qbasis[qrank][i]=bw[i]/n;
    qrank++;
  }
  std::array<double,TP_NE> left{};for(int i=0;i<TP_NE;i++)left[i]=f[i].pmid;
  for(int h=0;h<qrank;h++){double d=0;for(int i=0;i<TP_NE;i++)d+=left[i]*qbasis[h][i];for(int i=0;i<TP_NE;i++)left[i]-=d*qbasis[h][i];}
  double leftnorm=0;for(double y:left)leftnorm+=y*y;leftnorm=std::sqrt(leftnorm);
  if(leftnorm>1e-14){
    for(double&y:left)y/=leftnorm;
    I g;
    g=exact_taylor_combination(left,m,dx);
    /*
    double largest=-1;
    for(int j=0;j<TP_NV;j++){
      I gj;for(int i=0;i<TP_NE;i++)gj=gj+I(left[i])*f[i].grad[j];
      I contribution=gj*dx[j];g=g+contribution;
      if(contribution.width()>largest){largest=contribution.width();split_hint=j;}
    }
    */
    // Choose the widest coordinate; exact Taylor collection has already
    // cancelled all like linear and higher-order monomials.
    split_hint=0;for(int j=1;j<TP_NV;j++)if(b.x[j].width()>b.x[split_hint].width())split_hint=j;
    if(!g.zero()){lin++;return DROP;}
  }
  auto rows=select_rows(f);auto cols=select_cols(f,rows);
  double A[RK][RK]{},C[RK][RK]{};
  for(int i=0;i<RK;i++)for(int j=0;j<RK;j++)A[i][j]=f[rows[i]].gmid[cols[j]];
  if(!inverse(A,C)){invfail++;return KEEP;}
  std::array<bool,TP_NE> chosen{};for(int r:rows)chosen[r]=true;
  // Each omitted equation gives one left-null consistency test.
  for(int extra=0;extra<TP_NE;extra++)if(!chosen[extra]){
    double alpha[RK]{};
    for(int k=0;k<RK;k++)for(int j=0;j<RK;j++)alpha[k]+=f[extra].gmid[cols[j]]*C[j][k];
    I g=f[extra].pcenter;
    for(int k=0;k<RK;k++)g=g-I(alpha[k])*f[rows[k]].pcenter;
    for(int j=0;j<TP_NV;j++){
      I gj=f[extra].grad[j];
      for(int k=0;k<RK;k++)gj=gj-I(alpha[k])*f[rows[k]].grad[j];
      g=g+gj*dx[j];
    }
    if(!g.zero()){lin++;return DROP;}
  }
  I PJ[RK][RK]{},rhs[RK]{};
  std::array<bool,TP_NV> chosen_col{};for(int c:cols)chosen_col[c]=true;
  for(int i=0;i<RK;i++){
    for(int k=0;k<RK;k++)rhs[i]=rhs[i]-I(C[i][k])*f[rows[k]].pcenter;
    for(int j=0;j<TP_NV;j++)if(!chosen_col[j])
      for(int k=0;k<RK;k++)rhs[i]=rhs[i]-I(C[i][k])*f[rows[k]].grad[j]*dx[j];
    for(int j=0;j<RK;j++)for(int k=0;k<RK;k++)PJ[i][j]=PJ[i][j]+I(C[i][k])*f[rows[k]].grad[cols[j]];
  }
  std::array<I,RK>d{};for(int i=0;i<RK;i++)d[i]=dx[cols[i]];
  for(int pass=0;pass<2;pass++)for(int i=0;i<RK;i++){
    if(PJ[i][i].zero())continue;
    I rr=rhs[i];
    for(int j=0;j<RK;j++)if(i!=j)rr=rr-PJ[i][j]*d[j];
    d[i]=cap(d[i],rr/PJ[i][i]);if(d[i].empty()){nw++;return DROP;}
  }
  for(int i=0;i<RK;i++){int j=cols[i];b.x[j]=cap(b.x[j],I(m[j])+d[i]);if(b.x[j].empty()){nw++;return DROP;}}
  return KEEP;
}

}

int main(int argc,char**argv){
  std::uint64_t limit=1000000;int maxdepth=250,shard_index=0,shard_count=1;
  for(int i=1;i<argc;i++){
    std::string arg=argv[i];
    auto need=[&](){if(i+1>=argc){std::cerr<<"missing value after "<<arg<<"\n";std::exit(64);}return argv[++i];};
    if(arg=="--max-boxes")limit=std::strtoull(need(),nullptr,10);
    else if(arg=="--max-depth")maxdepth=std::atoi(need());
    else if(arg=="--shard-index")shard_index=std::atoi(need());
    else if(arg=="--shard-count")shard_count=std::atoi(need());
    else if(arg=="--help"){
      std::cout<<"usage: makeev_d5_interval [--max-boxes N] [--max-depth N] "
                 "[--shard-index I --shard-count N]\n"
                 "N=0 removes the box limit; shard count must be a power of two.\n";
      return 0;
    }else{std::cerr<<"unknown option: "<<arg<<"\n";return 64;}
  }
  if(limit==0)limit=std::numeric_limits<std::uint64_t>::max();
  if(shard_count<1||(shard_count&(shard_count-1))||shard_index<0||shard_index>=shard_count){
    std::cerr<<"invalid shard: require 0 <= index < a power-of-two count\n";return 64;
  }
  double B=1.0;Box first;for(auto&x:first.x)x=I(-B,B);
  for(int bit=0;(1<<bit)<shard_count;bit++){
    int variable=bit%TP_NV;double mid=first.x[variable].mid();
    if((shard_index>>bit)&1)first.x[variable].lo=mid;else first.x[variable].hi=mid;
  }
  std::vector<Box> stack{first};std::uint64_t boxes=0,range=0,lin=0,invfail=0,nw=0,unresolved=0;int deep=0;
  while(!stack.empty()&&boxes<limit){Box b=stack.back();stack.pop_back();boxes++;deep=std::max(deep,b.depth);
    int hint=-1;if(process(b,range,lin,invfail,nw,hint)==DROP)continue;
    if(b.depth>=maxdepth){unresolved++;continue;}
    int k=hint>=0?hint:0;
    if(hint<0)for(int j=1;j<TP_NV;j++)if(b.x[j].width()>b.x[k].width())k=j;
    double m=b.x[k].mid();Box l=b,r=b;l.depth=r.depth=b.depth+1;l.x[k].hi=m;r.x[k].lo=m;stack.push_back(r);stack.push_back(l);
  }
  std::cout<<"boxes "<<boxes<<" pending "<<stack.size()<<" range "<<range<<" consistency "<<lin
           <<" inverse_fail "<<invfail<<" newton "<<nw<<" unresolved "<<unresolved<<" deepest "<<deep<<"\n";
  if(stack.empty()&&unresolved==0){
    std::cout<<"Certificate verified for shard "<<shard_index<<"/"<<shard_count<<".\n";
    if(shard_count==1)std::cout<<"Certificate verified: no real tensor zero exists.\n";
    return 0;
  }
  std::cout<<"INCOMPLETE: increase --max-boxes or --max-depth.\n";
  return 2;
}
