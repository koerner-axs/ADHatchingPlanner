#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("avx2")

// #include <bits/stdc++.h>
// #include <stdio.h>
// #include <ext/pb_ds/assoc_container.hpp>
// #include <ext/pb_ds/tree_policy.hpp>

#include <vector>
#include <fstream>
#include <iostream>
#include <crtdbg.h>
#include <math.h>
#include <string>
#include <set>
//#include <stdio.h>

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


// Custom includes
#include "TinyPngOut-cpp/TinyPngOut.hpp"


const int MAXSIZE = 4096;

// temp
vector<vector<pii>> offsets;

/*
int tocover[maxsize][maxsize];
int canfarjump[maxsize][maxsize];



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

pii sample(int minjumpdist, int maxjumpdist) {
    int dist = rand() % (maxjumpdist - minjumpdist + 1);
    //int dist = minjumpdist;
    //cout << offsets[dist].size() << " " << ((dist+minjumpdist)*8) << endl;
    return {dist, rand() % ((dist+minjumpdist)*8)};
}

bool insiderect(pii point, pair<pii, pii> rect) {
    return point.fi >= rect.fi.fi && point.fi <= rect.se.fi && point.se >= rect.fi.se && point.se <= rect.se.se;
}

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






void save_matrix_as_png(vector<vector<int>> matrix, pair<int, int> size, string filename) {
    uint8_t* array = new uint8_t[size.fi*size.se*3];
    for (int x = 0; x < size.fi; x++) {
        for (int y = 0; y < size.se; y++) {
            uint8_t t = (uint8_t)matrix[x][y] * 0xFF;
            for (int i = 0; i < 3; i++) {
                array[x*size.fi*3+y*3+i] = t;
            }
        }
    }

    try {
		ofstream file;
        file.open(filename, ios::binary);
		TinyPngOut pngout((uint32_t)size.fi, (uint32_t)size.se, file);
		pngout.write(array, (size_t)(size.fi*size.se));
	} catch (const char *msg) {
		cout << "Failed to save matrix to file " << filename << endl;
        exit(1);
	}
}


bool isConflictFree(pair<int, int>& point, vector<pair<int, int>>& carry_forward,
                    vector<pair<int, int>>& cooling_down_points, int cooldown, const int mindist) {
    int a = point.fi;
    int b = point.se;
    for (int idx = (int)cooling_down_points.size()-1; idx >= 0 && cooldown > 0; idx--, --cooldown) {
        int c = cooling_down_points[idx].fi;
        int d = cooling_down_points[idx].se;
        if ((a-c)*(a-c) + (b-d)*(b-d) < mindist * mindist)
            return false;
    }
    for (int idx = (int)carry_forward.size()-1; idx >= 0 && cooldown > 0; idx--, --cooldown) {
        int c = carry_forward[idx].fi;
        int d = carry_forward[idx].se;
        if ((a-c)*(a-c) + (b-d)*(b-d) < mindist * mindist)
            return false;
    }
    return true;
}


struct TileEntry {
    int borderness;
    int posX, posY;
public:
    static bool comparison_func(TileEntry& lhs, TileEntry& rhs) {
        if (lhs.borderness != rhs.borderness) {
            return lhs.borderness < rhs.borderness;
        }
        if (lhs.posX != rhs.posX) {
            return lhs.posX < rhs.posX;
        }
        return lhs.posY < rhs.posY;
    }
};


class SpatialPartition {
private:
    vector<vector<set<TileEntry, TileEntry::comparison_func, vector<TileEntry>>>> data;
    pair<int, int> size;
    int tile_size;
public:
    SpatialPartition(pair<int, int> board_size, int tile_size);
    void build_partition(vector<vector<int>>* cell_state);
    int sample_in_range(pair<int, int>& out, int mindist, int maxdist, int maxiter);
    void update(pair<int, int> cell, bool remove);
};


struct Statistics {
    int initial_num_targets, initial_num_targets_longjump, initial_num_targets_nonlongjump;
    int current_num_targets, current_num_targets_longjump, current_num_targets_nonlongjump;
    int phase_one_num_tranches, phase_one_num_tranche_disappointees, phase_one_num_tranches_in_avg;
    double phase_one_avg_tranche_size;
};


class Board {
private:
    vector<vector<int>> cell_state;
    vector<vector<int>> cell_allow_farjump;
    pair<int, int> size;
    int farjump_threshold;
    int current_phase;
public:
    int minjumpdist, maxjumpdist, farjumpdepth, cooldown;
    vector<pair<int, int>> solution;
    Statistics statistics;
    SpatialPartition* partition;
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
    Board(istream& input);
    void initialize(int minjump, int maxjump, int farjump_borderdist, int cooldown_time);
    void solve();
    void print_solution_path();
    void print_board();
    void print_statistics();
    void save_board(string filename);
    static Board* read_board_from_file(string filename);
};


Board::Board(istream& input) {
    cell_state = vector<vector<int>>(MAXSIZE, vector<int>(MAXSIZE, 0));
    cell_allow_farjump = vector<vector<int>>(MAXSIZE, vector<int>(MAXSIZE, 0));

    // Read board size
    input >> size.fi >> size.se;

    _init_statistics();

    // Read board
    statistics.initial_num_targets = 0;
    for (int i = 1; i <= size.fi; i++) {
        for (int j = 1; j <= size.se; j++) {
            char c;
            input >> c;
            _ASSERT(c == '0' || c == '1');
            if (c == '1') {
                cell_state[i][j] = 1;
                statistics.initial_num_targets++;
            }
        }
    }
    statistics.current_num_targets = statistics.initial_num_targets;
}


void Board::initialize(int minjump, int maxjump, int farjump_borderdist, int cooldown_time) {
    minjumpdist = minjump;
    maxjumpdist = maxjump;
    farjumpdepth = farjump_borderdist;
    cooldown = cooldown_time;

    // Calculate the number of target cells in a neighborhood of given size for all cells.
    // If all cells in a neighborhood are target cells, we assume the center cell to be an
    // acceptable target for a long jump.
    _build_farjump();
}


void Board::solve() {
    solution.clear();
    current_phase = 1;


    save_board("./dbg/p0.png");


    // Phase one
    vector<pair<int, int>> carry_forward; // List of cells that need to cool down first
    int num_tranches = 0;
    while (current_phase == 1) {
        cout << "Tranche number " << (++num_tranches) << endl;
        print_statistics();

        vector<pair<int, int>> tranche = _build_tranche_phase_one(carry_forward);
        solution.insert(solution.end(), tranche.begin(), tranche.end());


        string name = string("./dbg/p1tranche") + to_string(num_tranches) + ".png";
        save_board(name);

        
        // Remember the last 'cooldown' many cells
        carry_forward.clear();
        for (int idx = max(0, (int)solution.size()-cooldown); idx < (int)solution.size(); idx++) {
            carry_forward.push_back(solution[idx]);
        }
    }

    // Phase two

    //cout << ((int)solution.size()) << endl;
    //exit(0);
    return;
}


void Board::print_solution_path() {
    cout << "Path:" << endl;
    for (int i = 0; i < (int)solution.size(); i++) {
        cout << solution[i].fi << " " << solution[i].se << endl;
    }
}


void Board::print_board() {
    cout << "Remaining:" << endl;
    for (int i = 1; i <= size.fi; i++) {
        for (int j = 1; j <= size.se; j++) {
            cout << cell_state[i][j] << " ";
            //cout << (tocover[i][j] + canfarjump[i][j]) << " ";
        }
        cout << endl;
    }
}

/*
int initial_num_targets, initial_num_targets_longjump, initial_num_targets_nonlongjump;
    int current_num_targets, current_num_targets_longjump, current_num_targets_nonlongjump;
    int phase_one_num_tranches, phase_one_num_tranche_disappointees, phase_one_num_tranches_in_avg;
    double phase_one_avg_tranche_size;*/

#define PRINT_NUM_AND_PERCENTOF(a,b) a << " (" << /*setprecision(2) <<*/ (100 * (a / (double)b)) << "%)"
void Board::print_statistics() {
    int area = size.fi * size.se;
    cout << "Current statistics:" << endl;
    cout << "Size of board: (" << size.fi << " x " << size.se << ")" << endl;
    cout << "Of which " << PRINT_NUM_AND_PERCENTOF(statistics.current_num_targets, area) << " cells are targets" << endl;
    cout << "   of which " << PRINT_NUM_AND_PERCENTOF(statistics.current_num_targets_longjump, statistics.current_num_targets) << " are eligible for a long jump, whereas " << PRINT_NUM_AND_PERCENTOF(statistics.current_num_targets_nonlongjump, statistics.current_num_targets) << " are not" << endl;
    cout << "to be extended" << endl;
}


void Board::save_board(string filename) {
    save_matrix_as_png(cell_state, size, filename);
}


Board* Board::read_board_from_file(string filename) {
    Board* board = nullptr;
    ifstream file;
    file.open(filename, ios::in);
    if (file.is_open()) {
        board = new Board(file);
        file.close();
    } else {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    return board;
}


void Board::_build_farjump() {
    statistics.initial_num_targets_longjump = 0;

    farjump_threshold = (2*farjumpdepth+1)*(2*farjumpdepth+1);
    for (int y = farjumpdepth; y < MAXSIZE - farjumpdepth; y++) {
        vi colcount(MAXSIZE, 0);
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
        for (int x = farjumpdepth+1; x < MAXSIZE - farjumpdepth; x++) {
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
    current_phase = phase_to_end + 1;
    //cout << "Exiting..." << endl;
    //exit(0);
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
        const int max_allowed_disappoint = 10;
        if (statistics.phase_one_num_tranche_disappointees > max_allowed_disappoint) {
            _declare_end_of_phase(1);
        }
    }

    if (!did_disappoint) { 
        statistics.phase_one_num_tranche_disappointees = max(0, statistics.phase_one_num_tranche_disappointees - 1);
        statistics.phase_one_avg_tranche_size = statistics.phase_one_avg_tranche_size * statistics.phase_one_num_tranches_in_avg + (double)tranche.size();
        statistics.phase_one_avg_tranche_size /= (double)(statistics.phase_one_num_tranches_in_avg + 1);
        statistics.phase_one_num_tranches_in_avg++;
    }

    // Update longjump and nonlongjump statistics
    for (pair<int, int> p : tranche) {
        statistics.current_num_targets--;
        if (cell_allow_farjump[p.fi][p.se] >= farjump_threshold) {
            statistics.current_num_targets_longjump--;
        } else {
            statistics.current_num_targets_nonlongjump--;
        }
    }

    /*set<pair<int, int>> s;
    pair<int, int> x;
    for (pair<int, int> p : tranche) {
        if (s.count(p)) {
            cout << "Point (" << p.fi << " " << p.se << ") was doubled" << endl;
            x = p;
        } else
        s.insert(p);
    }
    for (int idx = 0; idx < (int)tranche.size(); idx++) {
        if (tranche[idx] == x) {
            cout << "Doubling occured at index " << idx << endl;
            break;
        }
    }

    int lj = 0;
    int nlj = 0;
    rep(i,1,size.fi+1)
        rep(j,1,size.se+1) {
            if (cell_state[i][j]) {
                if (cell_allow_farjump[i][j] >= farjump_threshold)
                    lj++;
                else
                    nlj++;
            }
        }
    if (statistics.current_num_targets_longjump != lj || statistics.current_num_targets_nonlongjump != nlj) {
        print_statistics();
        cout << "Length of tranche " << ((int)tranche.size()) << endl;
        cout << "Longjump targets mismatched " << statistics.current_num_targets_longjump << " compared to " << lj << endl;
        cout << "Nonlongjump targets mismatched " << statistics.current_num_targets_nonlongjump << " compared to " << nlj << endl;
        exit(0);
    }*/
}


long long Board::_compute_max_time_phase_one() {
    //long long base_time_per_cell = 2 * ((2*maxjumpdist+1)*(2*maxjumpdist+1) - (2*minjumpdist+1)*(2*minjumpdist+1));
    long long base_time_per_cell = 4 * (8*minjumpdist - 4);
    double proportion_border_covered = 1.0 - (statistics.current_num_targets_nonlongjump / (double)statistics.initial_num_targets_nonlongjump);
    
    //cout << base_time_per_cell << " " << proportion_border_covered << endl;

    /*if (statistics.phase_one_num_tranches_in_avg >= 5) {
        return min(statistics.phase_one_avg_tranche_size, (double)statistics.current_num_targets_nonlongjump) * base_time_per_cell;
    } else {
        return statistics.current_num_targets_nonlongjump * base_time_per_cell;
    }*/

    
    if (statistics.phase_one_num_tranches_in_avg >= 5) {
        if (proportion_border_covered < 0.7) {
            return min(statistics.phase_one_avg_tranche_size, (double)statistics.current_num_targets_nonlongjump) * base_time_per_cell;
        } else if (proportion_border_covered < 0.9) {
            return min(statistics.phase_one_avg_tranche_size, (double)statistics.current_num_targets_nonlongjump) * base_time_per_cell * 4;
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
    
    pii currentpos = {-1, -1};
    /*for (int i = 1; i <= size.fi; i++) {
        for (int j = 1; j <= size.se; j++) {
            if (cell_state[i][j]) {
                currentpos = {i,j};
                break;
            }
        }
    }*/
    int row = rand() % size.fi;
    int col = rand() % size.se;
    while (currentpos.fi == -1) {
        for (int i = col; i < size.se; i++) {
            if (cell_state[row+1][i+1]) {
                currentpos = {row+1, i+1};
                break;
            }
        }
        if (currentpos.fi == -1) {
            for (int i = 0; i < col; i++) {
                if (cell_state[row+1][i+1]) {
                    currentpos = {row+1, i+1};
                    break;
                }
            }
        }
        row = (row + 1) % size.fi;
    }
    cell_state[currentpos.fi][currentpos.se] = 0;
    tranche.push_back(currentpos);

    //cout << max_time_spent << endl;


    const double failure_max_time_proportion = 0.001;
    const double resetfactor = 0.99;
    long long first_failure = 0;
    while (calc_time < max_time_spent) {



        pii s = sample(minjumpdist, maxjumpdist);
        s = offsets[s.fi][s.se];
        pii probe = {currentpos.fi + s.fi, currentpos.se + s.se};
        calc_time += 3;
        if (insiderect(probe, {{1, 1}, {size.fi, size.se}}) && cell_state[probe.fi][probe.se]) {
            calc_time += 8;
            /*if (cell_allow_farjump[probe.fi][probe.se] >= farjump_threshold) {
                if (rand() % 128 < 96) {
                    goto reject;
                }
            }*/
            if (isConflictFree(probe, carry_forward, tranche, cooldown, minjumpdist)) {
                tranche.pb(probe);
                //last_sites.pb(probe);
                //mark_keepout(probe, iteration);
                cell_state[probe.fi][probe.se] = 0;
                currentpos = probe;
                first_failure = calc_time; // Reset watchdog timer
                //if ((int)tranche.size() % 1000 == 0)
                //    cout << ((int)tranche.size()) << endl;
                //cout << currentpos.fi << " " << currentpos.se << endl;
            }
        }

        reject:

        // If for a long enough time no step could be made, reset by
        if ((calc_time-first_failure) >= (long long)ceil(failure_max_time_proportion * max_time_spent)) {
            int rollbackits = ceil(((int)tranche.size()) * (1-resetfactor));
            for (int i = 0; i < rollbackits; i++) {
                cell_state[tranche.back().fi][tranche.back().se] = 1;
                tranche.pop_back();
            }
            currentpos = tranche.back();
            first_failure = calc_time; // Reset watchdog timer
        }



    }

    _update_stats_phase_one_tranche(tranche);
    return tranche;
}


const int minjumpdist = 20;
const int maxjumpdist = 300;
//const int minjumpdist = 4;
//const int maxjumpdist = 50;
const int farjumpdepth = 3;
const int timeout = 10;

int main() {
    FIO;


    // temp
    build_offsets(minjumpdist, maxjumpdist);


    Board* board = Board::read_board_from_file("./inputs/large.txt");
    board->initialize(minjumpdist, maxjumpdist, farjumpdepth, timeout);
    board->solve();

    board->save_board("./last_output.png");
    //board->print_solution_path();
    //board->print_board();

    return 0;
}


/*
int main() {
    FIO;

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

    //cout << countones << " " << countzeros << " " << fillratio << endl;

    // rep(i,0,2) {
    //     rep(j,0,offsets[i].size()) {
    //         cout << "(" << offsets[i][j].fi << "," << offsets[i][j].se << ") ";
    //     }
    //     cout << endl;
    // }
    // return 0;

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

    return 0;
}
*/