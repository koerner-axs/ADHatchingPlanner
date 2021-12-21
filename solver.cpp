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
const int maxjumpdist = 300;
const int minjumpdist = 20;
const int timeout = 5;
const int farjumpdepth = 3;

/*
int tocover[maxsize][maxsize];
int canfarjump[maxsize][maxsize];
vector<vector<pii>> offsets;

void build_offsets(int min, int max) {
    rep(i,min,max+1) {
        vector<pii> vec;
        rep(j,-i,i+1) vec.pb({-i, j}), vec.pb({i, j});
        rep(j,-i+1,i) vec.pb({j, -i}), vec.pb({j, i});
        offsets.pb(vec);
        //cout << vec.size() << " ";
    }
    //cout << endl;
}

void build_canfarjump() {
    const int thresh = (2*farjumpdepth+1)*(2*farjumpdepth+1);
    for (int y = farjumpdepth; y < maxsize - farjumpdepth; y++) {
        vi colcount(maxsize, 0);
        int var = 0;
        for (int i = 0; i <= farjumpdepth*2; i++) {
            for (int j = -farjumpdepth; j <= farjumpdepth; j++) {
                colcount[i] += tocover[i][y+j];
            }
            var += colcount[i];
        }
        canfarjump[farjumpdepth][y] = var;// >= thresh;
        for (int x = farjumpdepth+1; x < maxsize - farjumpdepth; x++) {
            for (int j = -farjumpdepth; j <= farjumpdepth; j++) {
                colcount[x+farjumpdepth] += tocover[x+farjumpdepth][y+j];
            }
            var += colcount[x+farjumpdepth] - colcount[x-farjumpdepth-1];
            canfarjump[x][y] = var;// >= thresh;
        }
    }
}

pii sample() {
    //int dist = rand() % (maxjumpdist - minjumpdist + 1);
    int dist = minjumpdist;
    //cout << offsets[dist].size() << " " << ((dist+minjumpdist)*8) << endl;
    return {dist, rand() % ((dist+minjumpdist)*8)};
}

bool insiderect(pii point, pair<pii, pii> rect) {
    return point.fi >= rect.fi.fi && point.fi <= rect.se.fi && point.se >= rect.fi.se && point.se <= rect.se.se;
}

bool isallowed(pii point, int iteration, vector<pii>& steps) {
    for (int offset = 1; offset <= timeout; offset++) {
        int idx = (int)steps.size() - offset;
        if (idx < 0)
            break;
        pii p = steps[idx];
        //cout << idx << " " << p.fi << " " << p.se << endl;
        if (insiderect(point, {{p.fi - minjumpdist + 1, p.se - minjumpdist + 1}, {p.fi + minjumpdist - 1, p.se + minjumpdist - 1}}))
            return false;
    }
    return true;
}
*/

struct Statistics {
    int initial_num_targets, initial_num_targets_longjump, initial_num_targets_nonlongjump;
    int current_num_targets, current_num_targets_longjump, current_num_targets_nonlongjump;
    int phase_one_num_tranches, phase_one_num_tranche_disappointees, phase_one_num_tranches_in_avg;
    double phase_one_avg_tranche_size;
};

class Board {
private:
    int cell_state[maxsize][maxsize];
    int cell_allow_farjump[maxsize][maxsize];
    pair<int, int> size;
    int farjump_threshold;
public:
    int minjumpdist, maxjumpdist, farjumpdepth;
    Statistics statistics;
private:
    void _build_farjump();
    void _init_statistics();
    void _declare_end_of_phase(int phase_to_end);
    void _update_stats_phase_one_tranche(vector<pair<int, int>>& tranche);
    void _update_stats_phase_two_tranche(vector<pair<int, int>>& tranche);
    long long _compute_max_time_phase_one();
    long long _compute_max_time_phase_two();
    vector<pair<int, int>> _build_tranche_phase_one(vector<pair<int, int>>& carry_forward);
    vector<pair<int, int>> _build_tranche_phase_two(vector<pair<int, int>>& carry_forward);
public:
    Board(int minjumpdist, int maxjumpdist, int farjumpdepth);
    vector<pair<int, int>> solve();
};


Board::Board(int min, int max, int fjd) : minjumpdist(min), maxjumpdist(max), farjumpdepth(fjd) {
    // Read board size
    cin >> size.fi >> size.se;

    _init_statistics();

    // Read board
    statistics.initial_num_targets = 0;
    for (int i = 1; i <= size.fi; i++) {
        for (int j = 1; j <= size.se; j++) {
            char c;
            cin >> c;
            if (c == '1') {
                cell_state[i][j] = 1;
                statistics.initial_num_targets++;
            }
        }
    }
    statistics.current_num_targets = statistics.initial_num_targets;

    // Calculate the number of target cells in a neighborhood of given size for all cells.
    // If all cells in a neighborhood are target cells, we assume the center cell to be an
    // acceptable target for a long jump.
    _build_farjump();
}


vector<pair<int, int>> Board::solve() {
    return {};
}


void Board::_build_farjump() {
    statistics.initial_num_targets_longjump = 0;

    farjump_threshold = (2*farjumpdepth+1)*(2*farjumpdepth+1);
    for (int y = farjumpdepth; y < maxsize - farjumpdepth; y++) {
        vi colcount(maxsize, 0);
        int var = 0;
        for (int i = 0; i <= farjumpdepth*2; i++) {
            for (int j = -farjumpdepth; j <= farjumpdepth; j++) {
                colcount[i] += cell_state[i][y+j];
            }
            var += colcount[i];
        }
        cell_allow_farjump[farjumpdepth][y] = var;
        if (var >= farjump_threshold)
            statistics.initial_num_targets_longjump++;
        for (int x = farjumpdepth+1; x < maxsize - farjumpdepth; x++) {
            for (int j = -farjumpdepth; j <= farjumpdepth; j++) {
                colcount[x+farjumpdepth] += cell_state[x+farjumpdepth][y+j];
            }
            var += colcount[x+farjumpdepth] - colcount[x-farjumpdepth-1];
            cell_allow_farjump[x][y] = var;
            if (var >= farjump_threshold)
                statistics.initial_num_targets_longjump++;
        }
    }

    statistics.initial_num_targets_nonlongjump = statistics.initial_num_targets - statistics.initial_num_targets_longjump;
    statistics.current_num_targets_longjump = statistics.initial_num_targets_longjump;
    statistics.current_num_targets_nonlongjump = statistics.initial_num_targets_nonlongjump;
}


void Board::_init_statistics() {
    statistics.initial_num_targets = statistics.initial_num_targets_longjump = statistics.initial_num_targets_nonlongjump = 0;
    statistics.current_num_targets = statistics.current_num_targets_longjump = statistics.current_num_targets_nonlongjump = 0;
    statistics.phase_one_num_tranches = statistics.phase_one_num_tranche_disappointees = 0;
    statistics.phase_one_num_tranches_in_avg = 0;
    statistics.phase_one_avg_tranche_size = 0.0;
}


void Board::_declare_end_of_phase(int phase_to_end) {
    cout << "End of phase " << phase_to_end << endl;
    cout << "Exiting..." << endl;
    exit(0);
}


void Board::_update_stats_phase_one_tranche(vector<pair<int, int>>& tranche) {
    statistics.phase_one_num_tranches++;

    // Check if new tranche is of underwhelming size compared to running average. If so add one to number of consecutive disappointing
    // tranches. If that count exceeds a threshold, phase one is declared to have ended.
    bool did_disappoint = false;
    if (statistics.phase_one_num_tranches_in_avg > 5) { // Need to have a min number of components so that average is stable enough
        const double min_no_disappoint = 0.4; // Constant depends on the variability to be expected
        if ((double)tranche.size() < min_no_disappoint * statistics.phase_one_avg_tranche_size) {
            did_disappoint = true;
            statistics.phase_one_num_tranche_disappointees++;
        }
        const int max_allowed_disappoint = 5;
        if (statistics.phase_one_num_tranche_disappointees > max_allowed_disappoint) {
            _declare_end_of_phase(1);
        }
    }

    if (!did_disappoint) { 
        statistics.phase_one_avg_tranche_size = statistics.phase_one_avg_tranche_size * statistics.phase_one_num_tranches_in_avg + (double)tranche.size();
        statistics.phase_one_avg_tranche_size /= (double)(statistics.phase_one_num_tranches_in_avg + 1);
        statistics.phase_one_num_tranches_in_avg++;
    }

    // Update longjump and nonlongjump statistics
    for (pair<int, int> p : tranche) {
        if (cell_allow_farjump[p.fi][p.se] >= farjump_threshold) {
            statistics.current_num_targets_longjump--;
        } else {
            statistics.current_num_targets_nonlongjump--;
        }
    }
}


long long Board::_compute_max_time_phase_one() {
    long long base_time_per_cell = 2 * ((2*maxjumpdist+1)*(2*maxjumpdist+1) - (2*minjumpdist+1)*(2*minjumpdist+1));
    double proportion_border_covered = statistics.current_num_targets_nonlongjump / (double)statistics.initial_num_targets_nonlongjump;
    
    if (statistics.phase_one_num_tranches_in_avg >= 5) {
        if (proportion_border_covered < 0.7) {
            return statistics.phase_one_avg_tranche_size * base_time_per_cell;
        } else if (proportion_border_covered < 0.9) {
            return statistics.phase_one_avg_tranche_size * base_time_per_cell * 4;
        } else {
            return statistics.current_num_targets_nonlongjump * base_time_per_cell * 8;
        }
    } else {
        if (proportion_border_covered < 0.9) {
            return statistics.current_num_targets_nonlongjump * base_time_per_cell;
        } else {
            return statistics.current_num_targets_nonlongjump * base_time_per_cell * 8;
        }
    }
}


vector<pair<int, int>> Board::_build_tranche_phase_one(vector<pair<int, int>>& carry_forward) {
    vector<pair<int, int>> tranche;
    long long calc_time = 0;
    const long long max_time_spent = _compute_max_time_phase_one();
    while (calc_time < max_time_spent) {
        /* code */
    }
    _update_stats_phase_one_tranche(tranche);
    return tranche;
}


int main() {
    FIO;

    /*
    int sx, sy;
    cin >> sx >> sy;

    int countones = 0;
    for (int i = 1; i <= sx; i++) {
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
    build_canfarjump();
    */

    //cout << countones << " " << countzeros << " " << fillratio << endl;

    // rep(i,0,2) {
    //     rep(j,0,offsets[i].size()) {
    //         cout << "(" << offsets[i][j].fi << "," << offsets[i][j].se << ") ";
    //     }
    //     cout << endl;
    // }
    // return 0;

    /*
    deque<pii> last_sites;
    vector<vector<pii>> tranches;
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
    while (iteration <= 200000) {//countones*0.85) {//countones) {
        pii s = sample();
        s = offsets[s.fi][s.se];
        pii probe = {currentpos.fi + s.fi, currentpos.se + s.se};
        tries++;
        if (insiderect(probe, {{1, 1}, {sx, sy}}) && tocover[probe.fi][probe.se]) {
            tries += 8;
            if (isallowed(probe, iteration, steps)) {
                steps.pb(probe);
                //last_sites.pb(probe);
                //mark_keepout(probe, iteration);
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
    */

    /*
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
            //cout << (tocover[i][j] + canfarjump[i][j]) << " ";
        }
        cout << endl;
    }

    //vector<cluster> clusters;
    */

    return 0;
}