#include <vector>
#include <set>
#include <iostream>
#include <string>


struct TileEntry {
    int borderness;
    int posX, posY;
public:
    static bool comparison_func(const TileEntry& lhs, const TileEntry& rhs) {
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
    std::vector<std::vector<std::set<TileEntry, decltype(TileEntry::comparison_func)*>>> data;
    std::pair<int, int> size;
    int tile_size;
public:
    SpatialPartition(std::pair<int, int> board_size, int tile_size);
    void build_partition(std::vector<std::vector<int>>* cell_state);
    int sample_in_range(std::pair<int, int>& out, int mindist, int maxdist, int maxiter);
    void update(std::pair<int, int> cell, bool remove);
};


struct Statistics {
    int initial_num_targets, initial_num_targets_longjump, initial_num_targets_nonlongjump;
    int current_num_targets, current_num_targets_longjump, current_num_targets_nonlongjump;
    int phase_one_num_tranches, phase_one_num_tranche_disappointees, phase_one_num_tranches_in_avg;
    double phase_one_avg_tranche_size;
};


class Board {
private:
    std::vector<std::vector<int>> cell_state;
    std::vector<std::vector<int>> cell_allow_farjump;
    std::pair<int, int> size;
    int farjump_threshold;
    int current_phase;
public:
    int minjumpdist, maxjumpdist, farjumpdepth, cooldown;
    std::vector<std::pair<int, int>> solution;
    Statistics statistics;
    SpatialPartition* partition;
private:
    void _build_farjump();
    void _init_statistics();
    void _declare_end_of_phase(int phase_to_end);
    void _update_stats_phase_one_tranche(std::vector<std::pair<int, int>>& tranche);
    void _update_stats_phase_two_tranche(std::vector<std::pair<int, int>>& tranche);
    long long _compute_max_time_phase_one();
    long long _compute_max_time_phase_two();
    std::vector<std::pair<int, int>> _build_tranche_phase_one(std::vector<std::pair<int, int>>& carry_forward);
    std::vector<std::pair<int, int>> _build_tranche_phase_two(std::vector<std::pair<int, int>>& carry_forward);
public:
    Board(std::istream& input);
    void initialize(int minjump, int maxjump, int farjump_borderdist, int cooldown_time);
    void solve();
    void print_solution_path();
    void print_board();
    void print_statistics();
    void save_board(std::string filename);
    static Board* read_board_from_file(std::string filename);
};
