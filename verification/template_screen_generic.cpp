#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>

using boost::multiprecision::cpp_int;
using Vec = std::vector<int>;
static int QB=13;

static long long gcdll(long long a,long long b){while(b){long long t=a%b;a=b;b=t;}return a;}
static long long lcmll(long long a,long long b){return a/gcdll(a,b)*b;}

struct Group {
    Vec inv;
    int m=1;
    std::vector<Vec> elems;
    std::map<Vec,int> eidx;
    std::vector<uint32_t> submasks;
    std::map<int,std::vector<int>> by_orbit;
    std::vector<std::vector<uint8_t>> trivial;
    std::vector<int> pairidx;
};

static void gen_elems_rec(const Vec& inv,int i,Vec& cur,std::vector<Vec>& out){
    if(i==(int)inv.size()){out.push_back(cur);return;}
    for(int x=0;x<inv[i];++x){cur[i]=x;gen_elems_rec(inv,i+1,cur,out);} }

static Vec addv(const Vec&a,const Vec&b,const Vec&inv){Vec c(a.size());for(size_t i=0;i<a.size();++i)c[i]=(a[i]+b[i])%inv[i];return c;}

static uint32_t span_mask(const Group& G,const std::vector<int>& gens){
    std::set<Vec> S; Vec z(G.inv.size(),0); S.insert(z);
    for(int gi:gens){
        std::vector<Vec> prior(S.begin(),S.end());
        Vec g=G.elems[gi];
        int ord=1; Vec t=g;
        while(t!=z){t=addv(t,g,G.inv);++ord;}
        for(auto &x:prior){Vec y=x;for(int c=1;c<ord;++c){y=addv(y,g,G.inv);S.insert(y);}}
    }
    uint32_t mask=0; for(auto &x:S) mask |= (uint32_t(1)<<G.eidx.at(x)); return mask;
}

static void comb_rec(int n,int k,int start,std::vector<int>&cur,const std::function<void(const std::vector<int>&)>&f){
    if((int)cur.size()==k){f(cur);return;}
    for(int i=start;i<n;++i){cur.push_back(i);comb_rec(n,k,i+1,cur,f);cur.pop_back();}
}

static Group make_group(const Vec& inv,bool legacy){
    Group G;G.inv=inv;for(int d:inv)G.m*=d;
    Vec cur(inv.size(),0);gen_elems_rec(inv,0,cur,G.elems);for(int i=0;i<G.m;++i)G.eidx[G.elems[i]]=i;
    std::set<uint32_t> subs; int maxgens=legacy?std::min(3,(int)inv.size()):(int)inv.size();
    for(int r=0;r<=maxgens;++r){std::vector<int> c;comb_rec(G.m,r,0,c,[&](const std::vector<int>& gs){subs.insert(span_mask(G,gs));});}
    G.submasks.assign(subs.begin(),subs.end());
    std::sort(G.submasks.begin(),G.submasks.end(),[](uint32_t a,uint32_t b){int ca=__builtin_popcount(a),cb=__builtin_popcount(b);return ca==cb?a<b:ca<cb;});
    for(int si=0;si<(int)G.submasks.size();++si){int sz=__builtin_popcount(G.submasks[si]);G.by_orbit[G.m/sz].push_back(si);}
    long long L=1;for(int d:inv)L=lcmll(L,d);
    G.trivial.assign(G.submasks.size(),std::vector<uint8_t>(G.m,0));
    for(int si=0;si<(int)G.submasks.size();++si){
        uint32_t sm=G.submasks[si];
        for(int ci=0;ci<G.m;++ci){bool ok=true;auto &c=G.elems[ci];
            for(int hi=0;hi<G.m && ok;++hi) if(sm&(uint32_t(1)<<hi)){
                auto &h=G.elems[hi]; long long s=0;for(size_t j=0;j<inv.size();++j)s += (long long)c[j]*h[j]*(L/inv[j]);
                if(((s%L)+L)%L!=0)ok=false;
            }
            G.trivial[si][ci]=ok;
        }
    }
    G.pairidx.resize(G.m);
    for(int ci=0;ci<G.m;++ci){Vec cc(inv.size());for(size_t j=0;j<inv.size();++j){int x=(-QB*G.elems[ci][j])%inv[j];if(x<0)x+=inv[j];cc[j]=x;}G.pairidx[ci]=G.eidx.at(cc);}
    return G;
}

static cpp_int ipow(cpp_int a,int e){cpp_int r=1;while(e){if(e&1)r*=a;a*=a;e>>=1;}return r;}
static cpp_int gauss_binom(int Q,int d,int t){cpp_int num=1,den=1;for(int i=0;i<t;++i){num*=ipow(Q,d-i)-1;den*=ipow(Q,i+1)-1;}return num/den;}
static cpp_int lagrangians(int q,int d){cpp_int r=1;for(int i=1;i<=d/2;++i)r*=ipow(q,2*i-1)+1;return r;}
static cpp_int paired_options(int Q,int d){cpp_int r=0;for(int t=0;t<=d;++t)r+=gauss_binom(Q,d,t);return r;}
static std::string cstr(const cpp_int&x){return x.convert_to<std::string>();}
static std::string vstr(const Vec&v){std::ostringstream o;o<<"(";for(size_t i=0;i<v.size();++i){if(i)o<<",";o<<v[i];}o<<")";return o.str();}

static void shapes_rec(const std::vector<int>&sizes,int i,int rem,Vec&cur,std::vector<Vec>&out){
    if(rem==0){out.push_back(cur);return;} if(i==(int)sizes.size())return; int s=sizes[i];
    for(int k=0;k*s<=rem;++k){for(int z=0;z<k;++z)cur.push_back(s);shapes_rec(sizes,i+1,rem-k*s,cur,out);for(int z=0;z<k;++z)cur.pop_back();}
}

struct ProfileKey {
    cpp_int fam; Vec inv,shape,sp,pr;
};
struct ProfileCmp {bool operator()(ProfileKey const&a,ProfileKey const&b) const{
    if(a.fam!=b.fam)return a.fam<b.fam; if(a.inv!=b.inv)return a.inv<b.inv; if(a.shape!=b.shape)return a.shape<b.shape;if(a.sp!=b.sp)return a.sp<b.sp;return a.pr<b.pr;}};

struct Results {uint64_t viable=0;std::map<ProfileKey,uint64_t,ProfileCmp> profiles;std::set<cpp_int> famsizes;};

static void enumerate_template(const Group&G,const Vec&shape,Results&R){
    std::map<int,int> counts;for(int s:shape)counts[s]++;
    std::vector<std::pair<std::vector<int>,int>> cats;for(auto &[s,k]:counts)cats.push_back({G.by_orbit.at(s),k});
    std::vector<int> mult(G.m,0);
    std::function<void(int)> catrec;
    catrec=[&](int ci){
        if(ci==(int)cats.size()){
            std::vector<char> done(G.m,0);Vec sp,pr;
            for(int c=0;c<G.m;++c)if(!done[c]){int cc=G.pairidx[c];if(mult[c]!=mult[cc])return;if(cc==c){if(mult[c]&1)return;if(mult[c])sp.push_back(mult[c]);done[c]=1;}else{if(mult[c])pr.push_back(mult[c]);done[c]=done[cc]=1;}}
            std::sort(sp.begin(),sp.end());std::sort(pr.begin(),pr.end());cpp_int fam=1;for(int d:sp)fam*=lagrangians(QB,d);for(int d:pr)fam*=paired_options(QB*QB,d);
            ProfileKey key{fam,G.inv,shape,sp,pr};R.viable++;R.profiles[key]++;R.famsizes.insert(fam);return;
        }
        auto &opts=cats[ci].first;int k=cats[ci].second;
        std::function<void(int,int)> choose=[&](int start,int left){
            if(left==0){catrec(ci+1);return;}
            for(int oi=start;oi<(int)opts.size();++oi){int si=opts[oi];for(int c=0;c<G.m;++c)mult[c]+=G.trivial[si][c];choose(oi,left-1);for(int c=0;c<G.m;++c)mult[c]-=G.trivial[si][c];}
        };
        choose(0,k);
    };
    catrec(0);
}

int main(int argc,char**argv){bool legacy=false;for(int i=1;i<argc;++i){std::string a=argv[i];if(a=="--legacy")legacy=true;else QB=std::stoi(a);}
    std::vector<Vec> allgroups={{2},{3},{2,2},{4},{5},{6},{7},{2,2,2},{2,4},{8},{3,3},{9},{10},{11},{2,6},{12},{13},{14},{15},{2,2,2,2},{2,2,4},{2,8},{4,4},{16},{17},{3,6},{18}};
    std::vector<Vec> groups; for(auto g:allgroups){int m=1;for(int d:g)m*=d;if(m%QB!=0 && (QB*QB-1)%(*std::max_element(g.begin(),g.end()))==0)groups.push_back(g);}
    Results R;
    for(auto inv:groups){Group G=make_group(inv,legacy);std::vector<int> sizes;for(auto &[s,x]:G.by_orbit)sizes.push_back(s);std::sort(sizes.rbegin(),sizes.rend());Vec cur;std::vector<Vec> shapes;shapes_rec(sizes,0,18,cur,shapes);
        uint64_t before=R.viable;for(auto &sh:shapes)enumerate_template(G,sh,R);std::cerr<<vstr(inv)<<" subgroups="<<G.submasks.size()<<" shapes="<<shapes.size()<<" viable="<<(R.viable-before)<<"\n";
    }
    std::cout<<(legacy?"LEGACY":"CORRECTED")<<" q="<<QB<<"\n";
    std::cout<<"viable_templates="<<R.viable<<"\nprofiles="<<R.profiles.size()<<"\ndistinct_family_sizes="<<R.famsizes.size()<<"\n";
    int rank=0;cpp_int target = (QB==13 ? cpp_int("12252982784") : cpp_int("5159780352"));for(auto &f:R.famsizes){rank++;if(rank<=10)std::cout<<"family_rank_"<<rank<<"="<<f<<"\n";if(f==target)std::cout<<"searched_family_distinct_rank="<<rank<<"\n";}
    int shown=0;for(auto &[k,c]:R.profiles){if(shown++<15 || (k.inv==Vec{3,3} && k.shape==Vec{9,9})){
        std::cout<<"PROFILE family="<<k.fam<<" H="<<vstr(k.inv)<<" shape="<<vstr(k.shape)<<" sp="<<vstr(k.sp)<<" pr="<<vstr(k.pr)<<" count="<<c;
        if(k.inv==Vec{3,3}&&k.shape==Vec{9,9})std::cout<<" SEARCHED";std::cout<<"\n";}
    }
}
