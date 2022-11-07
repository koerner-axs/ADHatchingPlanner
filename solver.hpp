#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <chrono>


#define fi first
#define se second
#define cout std::cout
#define endl std::endl

#ifdef _DEBUG
#define _DBG_TIME_PROCEDURE_A(x) cout << "enter procedure " << x << endl; auto time_measurement_start = std::chrono::steady_clock::now();
#define _DBG_TIME_PROCEDURE_B(x) auto time_measurement_end = std::chrono::steady_clock::now(); cout << "leave procedure " << x << " after " << (std::chrono::duration_cast<std::chrono::milliseconds>(time_measurement_end - time_measurement_start).count()) << " ms" << endl;
#else
#define _DBG_TIME_PROCEDURE_A(x)
#define _DBG_TIME_PROCEDURE_B(x)
#endif


struct TileEntry {
    int borderness;
    int posX, posY;
public:
    TileEntry(int borderness, int posX, int posY) : borderness(borderness), posX(posX), posY(posY) {}
    static bool comparison_func(const TileEntry& lhs, const TileEntry& rhs) {
        /*if (lhs.borderness != rhs.borderness) {
            return lhs.borderness < rhs.borderness;
        }*/
        if (lhs.posX != rhs.posX) {
            return lhs.posX < rhs.posX;
        }
        return lhs.posY < rhs.posY;
    }
    std::pair<int, int> pos_pair() const {
        return {posX, posY};
    }
};

typedef std::set<TileEntry, decltype(TileEntry::comparison_func)*> Tile;

class Board;
class SpatialPartition {
public:
    std::vector<std::vector<Tile>> data;
    std::pair<int, int> board_size;
    std::pair<int, int> size;
    int stride;
    Board* board;
private:
    void _convert_to_internal(std::pair<int, int> cell, std::pair<int, int>& part_cell, std::pair<int, int>& cell_offset);
    long long _find_allowed_in_tile(Tile& tile, std::vector<std::pair<int, int>>& allowed_targets, const int min_poss_borderness, std::pair<int, int>& current_best_cell, int& current_best_borderness);
public:
    SpatialPartition(Board* board, int tile_size);
    long long sample_in_range(std::pair<int, int> cell, TileEntry& out, int mindist, int maxdist, long long maxiter);
    void insert(std::pair<int, int> point);
    void insert(TileEntry entry);
    void remove(std::pair<int, int> point);
    void remove(TileEntry entry);
};


struct Statistics {
    int initial_num_targets, initial_num_targets_longjump, initial_num_targets_nonlongjump;
    int current_num_targets, current_num_targets_longjump, current_num_targets_nonlongjump;
    int phase_one_num_tranches, phase_one_num_tranche_disappointees, phase_one_num_tranches_in_avg;
    double phase_one_avg_tranche_size;
};


class Board {
public:
    std::vector<std::vector<int>> cell_state;
    std::vector<std::vector<int>> cell_borderness;
    std::vector<std::pair<int, int>> current_tranche;
    int farjump_threshold;
    int current_phase;
    std::pair<int, int> size;
    int minjumpdist, maxjumpdist, farjumpdepth, cooldown;
    std::vector<std::pair<int, int>> solution;
    Statistics statistics;
    SpatialPartition* partition;
private:
    void _build_borderness();
    void _build_spatial_partition();
    void _init_statistics();
    void _declare_end_of_phase(int phase_to_end);
    void _update_stats_phase_one_tranche();
    void _update_stats_phase_two_tranche();
    long long _compute_max_time_phase_one();
    long long _compute_max_time_phase_two();
    void _build_tranche_phase_one();
    void _build_tranche_phase_two();
public:
    Board(std::istream& input);
    void initialize(int minjump, int maxjump, int farjump_borderdist, int cooldown_time);
    void solve();
    void print_solution_path();
    void print_board();
    void print_statistics();
    void save_board(std::string filename);
    bool query_conflict_free(const std::pair<int, int>& probe);
    static Board* read_board_from_file(std::string filename);
};
