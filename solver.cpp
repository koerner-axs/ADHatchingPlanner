#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("avx2")

#include <bits/stdc++.h>
// #include <stdio.h>
// #include <ext/pb_ds/assoc_container.hpp>
// #include <ext/pb_ds/tree_policy.hpp>

using namespace std;
// using namespace __gnu_pbds;

#define FIO ios_base::sync_with_stdio(false); cin.tie(0);
#define mp make_pair
#define pb push_back
#define fi first
#define se second
#define rep(i, a, b) for(int i = a; i < (b); i++)
#define sz(x) (int)(x).size()
#define all(x) (x).begin(), (x).end()

using ll = long long;
using ld = long double;
using str = string;
using pii = pair<int, int>;
using pll = pair<ll, ll>;
using vi = vector<int>;
using vvi = vector<vector<int>>;
// using oset = tree<int,null_type,less<int>,rb_tree_tag,tree_order_statistics_node_update>;

const int maxsize = 4096;
const int maxjumpdist = 75;
const int minjumpdist = 20;
const int timeout = 5;

int tocover[maxsize][maxsize];
int keepout[maxsize][maxsize];
vector<vector<pii>> offsets;

void build_offsets(int min, int max) {
    rep(i,min,max+1) {
        vector<pii> vec;
        rep(j,-i,i+1) vec.pb({-i, j}), vec.pb({i, j});
        rep(j,-i+1,i) vec.pb({j, -i}), vec.pb({j, i});
        offsets.pb(vec);
        cout << vec.size() << " ";
    }
    cout << endl;
}

pii sample() {
    int dist = rand() % (maxjumpdist - minjumpdist + 1);
    //cout << offsets[dist].size() << " " << ((dist+minjumpdist)*8) << endl;
    return {dist, rand() % ((dist+minjumpdist)*8)};
}

bool insiderect(pii point, pair<pii, pii> rect) {
    return point.fi >= rect.fi.fi && point.fi <= rect.se.fi && point.se >= rect.fi.se && point.se <= rect.se.se;
}

bool isallowed(pii point, int iteration) {
    if (keepout[point.fi][point.se] == 0 || iteration <= keepout[point.fi][point.se])
        return true;
    return iteration - keepout[point.fi][point.se] >= timeout;
}

void mark_keepout(pii point, int iteration) {
    int left = max(1, point.fi - minjumpdist + 1);
    int right = min(maxsize, point.fi + minjumpdist - 1);
    int top = max(1, point.se - minjumpdist + 1);
    int bottom = min(maxsize, point.se + minjumpdist - 1);
    rep(i,left,right+1) {
        rep(j,top,bottom+1) {
            keepout[i][j] = iteration;
        }
    }
}


int main() {
    FIO;

    int sx, sy;
    cin >> sx >> sy;

    int countones = 0;
    for (int i = 1; i <= sy; i++) {
        for (int j = 1; j <= sy; j++) {
            char c;
            cin >> c;
            if (c == '1') {
                tocover[i][j] = 1;
                countones++;
            }
        }
    }
    int countzeros = sx * sy - countones;
    double fillratio = countones / (double)(sx * sy);
    build_offsets(minjumpdist, maxjumpdist);

    cout << countones << " " << countzeros << " " << fillratio << endl;

    // rep(i,0,2) {
    //     rep(j,0,offsets[i].size()) {
    //         cout << "(" << offsets[i][j].fi << "," << offsets[i][j].se << ") ";
    //     }
    //     cout << endl;
    // }
    // return 0;


    deque<pii> last_sites;
    vector<pii> steps;
    pii currentpos;
    rep(i,0,sx) {
        rep(j,0,sy) {
            if (tocover[i][j]) {
                currentpos = {i+1,j+1};
                last_sites.pb(currentpos);
                break;
            }
        }
    }
    int iteration = 1;
    int tries = 0;
    int maxtries = ((2*maxjumpdist-1)*(2*maxjumpdist-1) - (2*minjumpdist-1)*(2*minjumpdist-1)) * 16;
    const double resetfactor = 0.99;
    while (iteration <= countones*0.25) {//countones) {
        pii s = sample();
        s = offsets[s.fi][s.se];
        pii probe = {currentpos.fi + s.fi, currentpos.se + s.se};
        tries++;
        if (insiderect(probe, {{1, 1}, {sx, sy}}) && tocover[probe.fi][probe.se]) {
            tries += 8;
            if (isallowed(probe, iteration)) {
                steps.pb(probe);
                //last_sites.pb(probe);
                mark_keepout(probe, iteration);
                tocover[probe.fi][probe.se] = 0;
                currentpos = probe;
                iteration++;
                tries = 0;
                cout << iteration << endl;
                //cout << currentpos.fi << " " << currentpos.se << endl;
            }
        }
        if (tries >= maxtries) {
            int rollbackits = ceil(iteration * (1-resetfactor));
            iteration -= rollbackits;
            rep(i,0,rollbackits) {
                tocover[steps.back().fi][steps.back().se] = 1;
                steps.pop_back();
            }
            currentpos = steps.back();
            tries = 0;
        }
    }

    cout << endl;
    cout << endl;
    cout << endl;

    cout << "Path:" << endl;
    rep(i,0,(int)steps.size()) {
        cout << steps[i].fi << " " << steps[i].se << endl;
    }

    cout << "Remaining:" << endl;
    rep(i,1,sx+1) {
        rep(j,1,sy+1) {
            cout << tocover[i][j] << " ";
        }
        cout << endl;
    }

    //vector<cluster> clusters;


    return 0;
}