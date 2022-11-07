#include "solver.hpp"

#include <math.h>
#include <algorithm>
#include <random>


SpatialPartition::SpatialPartition(Board* board, int tile_size)
    : board(board), board_size(board->size), stride(tile_size) {
    
    size = {(int)std::ceil(board_size.fi / (double)stride), (int)std::ceil(board_size.se / (double)stride)};
    data = std::vector<std::vector<Tile>>(size.fi, std::vector<Tile>(size.se, Tile(TileEntry::comparison_func)));
}


std::pair<int, int> _limit_to_rect(std::pair<int, int> a, const std::pair<int, int> tl, const std::pair<int, int> br) {
    return {std::min(br.fi, std::max(tl.fi, a.fi)), std::min(br.se, std::max(tl.se, a.se))};
}


std::random_device rd;
std::mt19937 rng(rd());


long long SpatialPartition::sample_in_range(std::pair<int, int> cell, TileEntry& out, int mindist, int maxdist, long long maxiter) {
    std::pair<int, int> center_tile, offset;
    _convert_to_internal(cell, center_tile, offset);

    std::pair<int, int> rect_tl = { 1, 1 };
    std::pair<int, int> rect_br = board_size;
    std::pair<int, int> outer_tl_tile, outer_br_tile;
    _convert_to_internal(_limit_to_rect({ cell.fi - maxdist, cell.se - maxdist }, rect_tl, rect_br), outer_tl_tile, offset);
    _convert_to_internal(_limit_to_rect({ cell.fi + maxdist, cell.se + maxdist }, rect_tl, rect_br), outer_br_tile, offset);

    std::pair<int, int> inner_tl_tile, inner_br_tile;
    _convert_to_internal(_limit_to_rect({ cell.fi - mindist, cell.se - mindist }, rect_tl, rect_br), inner_tl_tile, offset);
    _convert_to_internal(_limit_to_rect({ cell.fi + mindist, cell.se + mindist }, rect_tl, rect_br), inner_br_tile, offset);

    std::pair<int, int> current_best_cell = { -1, -1 };
    const int minimum_possible_borderness = 0;
    int current_best_borderness = INT_MAX;
    long long cost_cum = 0;
    std::vector<std::pair<int, int>> candidates;

    // Random number to choose index
    int width = outer_br_tile.first - outer_tl_tile.first + 1;
    int height = outer_br_tile.second - outer_tl_tile.second + 1;
    int num_simple_samples = std::min(10, width * height); // Unlikely error mitigation

    // Simple sampling:
    // Randomly select a tile, but check first whether it is in the covered set.
    // This is done for the first few tiles as it is very fast.
    std::set<std::pair<int, int>> covered;
    std::uniform_int_distribution x_dist(outer_tl_tile.first, outer_br_tile.first);
    std::uniform_int_distribution y_dist(outer_tl_tile.second, outer_tl_tile.second);
    int iteration = 0;
    while (iteration < num_simple_samples) {
        //int x = outer_tl_tile.first + (rand() % width);
        //int y = outer_tl_tile.second + (rand() % height);
        int x = x_dist(rng);
        int y = y_dist(rng);

        if (covered.count({ x, y }))
            continue;
        else
            iteration++;

        cost_cum += _find_allowed_in_tile(data[x][y], candidates, minimum_possible_borderness, current_best_cell, current_best_borderness);
        if (!candidates.empty())
            goto success;

        covered.insert({ x, y });
    }

    // Explicit random permutation based sampling
    {
        // Build index vector and shuffle it.
        std::vector<std::pair<int, int>> indices(width * height, { -1, -1 });
        int counter = 0;
        for (int x = outer_tl_tile.first; x <= outer_br_tile.first; x++) {
            for (int y = outer_tl_tile.second; y <= outer_br_tile.second; y++) {
                indices[counter++] = { x, y };
            }
        }
        std::shuffle(indices.begin(), indices.end(), rng);

        // Sample by iterating through the permuted vector.
        // Filter out already covered tiles
        for (std::pair<int, int> tile : indices) {
            if (covered.count(tile))
                continue;
            cost_cum += _find_allowed_in_tile(data[tile.first][tile.second], candidates, minimum_possible_borderness, current_best_cell, current_best_borderness);
            if (!candidates.empty())
                goto success;
        }
    }

    /*
    int x_start = outer_tl_tile.first + (rand() % width);
    int y_start = outer_tl_tile.second + (rand() % height);

    // From start to the end of the scan field.
    for (int x = x_start; x <= outer_br_tile.first; x++) {
        int y;
        if (x == x_start)
            // Only scan from (including) starting point.
            y = y_start;
        else
            // Go full width of the scan field.
            y = outer_tl_tile.second;
        for (; y <= outer_br_tile.second; y++) {
            cost_cum += _find_allowed_in_tile(data[x][y], candidates, minimum_possible_borderness, current_best_cell, current_best_borderness);
            if (!candidates.empty())
                goto success;
        }
    }

    // From beginning of the scan field to the start of the scan.
    for (int x = outer_tl_tile.first; x <= x_start; x++) {
        int y_end;
        if (x == x_start)
            // Only scan to the starting point.
            y_end = y_start - 1; // -1 due to for loop condition <=
        else
            // Go full width of the scan field.
            y_end = outer_br_tile.second;
        for (int y = outer_tl_tile.second; y <= y_end; y++) {
            cost_cum += _find_allowed_in_tile(data[x][y], candidates, minimum_possible_borderness, current_best_cell, current_best_borderness);
            if (!candidates.empty())
                goto success;
        }
    }*/

    /*for (int x = outer_tl_tile.first; x <= outer_br_tile.first; x++) {
        for (int y = outer_tl_tile.second; y <= outer_br_tile.second; y++) {
            cost_cum += _find_allowed_optimum_in_tile(data[x][y], candidates, minimum_possible_borderness, current_best_cell, current_best_borderness);
        }
    }*/

    if (candidates.empty()) {
        // Did not find a viable target.
        out = TileEntry(current_best_borderness, current_best_cell.first, current_best_cell.second);
        return cost_cum;
    }

success:
    current_best_cell = candidates[rand() % ((int)candidates.size())];
    current_best_borderness = board->cell_borderness[current_best_cell.first][current_best_cell.second];

    out = TileEntry(current_best_borderness, current_best_cell.first, current_best_cell.second);
    return cost_cum;

    /*
    // The scanning pattern is the following with:
    // 0 = uncovered    1 = left    2 = top    3 = bottom    4 = right
    //
    // 0 0 0 0 0 0 0 0
    // 0 1 1 2 2 4 4 0
    // 0 1 1 2 2 4 4 0
    // 0 1 1 0 0 4 4 0
    // 0 1 1 0 0 4 4 0
    // 0 1 1 3 3 4 4 0
    // 0 1 1 3 3 4 4 0
    // 0 0 0 0 0 0 0 0

    //cout << outer_tl_tile.first << " " << outer_tl_tile.second << endl;
    //cout << outer_br_tile.first << " " << outer_br_tile.second << endl;

    // If none found write nothing to out, out is already {-1, -1}
    // Top strip
    if (current_best_borderness > minimum_possible_borderness) {
        for (int x = outer_tl_tile.first; x <= outer_br_tile.first; x++) {
            for (int y = outer_tl_tile.second; y < inner_tl_tile.second; y++) {
                cost_cum += _find_allowed_optimum_in_tile(data[x][y], minimum_possible_borderness, current_best_cell, current_best_borderness);
            }
        }
    }
    // Left strip
    if (current_best_borderness > minimum_possible_borderness) {
        for (int x = outer_tl_tile.first; x < inner_tl_tile.first; x++) {
            for (int y = inner_tl_tile.second; y <= inner_br_tile.second; y++) {
                cost_cum += _find_allowed_optimum_in_tile(data[x][y], minimum_possible_borderness, current_best_cell, current_best_borderness);
            }
        }
    }
    // Right strip
    if (current_best_borderness > minimum_possible_borderness) {
        for (int x = outer_br_tile.first; x > inner_br_tile.first; x--) {
            for (int y = inner_tl_tile.second; y <= inner_br_tile.second; y++) {
                cost_cum += _find_allowed_optimum_in_tile(data[x][y], minimum_possible_borderness, current_best_cell, current_best_borderness);
            }
        }
    }
    // Bottom strip
    if (current_best_borderness > minimum_possible_borderness) {
        for (int x = outer_tl_tile.first; x <= outer_br_tile.first; x++) {
            for (int y = outer_br_tile.second; y > inner_br_tile.second; y--) {
                cost_cum += _find_allowed_optimum_in_tile(data[x][y], minimum_possible_borderness, current_best_cell, current_best_borderness);
            }
        }
    }

    //for (int x = std::max(0, outer_tl_tile.first - 1); x <= outer_br_tile.first + 1; x++) {
    //    for (int y = std::max(0, outer_tl_tile.second - 1); y <= outer_br_tile.second + 1; y++) {
    //        cout << " " << dbg[x][y];
    //    }
    //    cout << endl;
    //}
    //exit(0);
    */
}


void SpatialPartition::insert(std::pair<int, int> point) {
    insert(TileEntry(board->cell_borderness[point.first][point.second], point.first, point.second));
}


void SpatialPartition::insert(TileEntry entry) {
    std::pair<int, int> tile, offset;
    _convert_to_internal({entry.posX, entry.posY}, tile, offset);
    data[tile.fi][tile.se].insert(entry);
}


void SpatialPartition::remove(std::pair<int, int> point) {
    remove(TileEntry(board->cell_borderness[point.first][point.second], point.first, point.second));
}


void SpatialPartition::remove(TileEntry entry) {
    std::pair<int, int> tile, offset;
    _convert_to_internal({entry.posX, entry.posY}, tile, offset);
    data[tile.fi][tile.se].erase(entry);
}


void SpatialPartition::_convert_to_internal(std::pair<int, int> cell,
    std::pair<int, int>& part_cell, std::pair<int, int>& cell_offset) {

    part_cell = {(int)std::floor(cell.fi / (double)stride), (int)std::floor(cell.se / (double)stride)};
    cell_offset = {cell.fi - (part_cell.fi * stride), cell.se - (part_cell.se * stride)};
}

int num_calls = 0;

long long SpatialPartition::_find_allowed_in_tile(Tile& tile,
    std::vector<std::pair<int, int>>& allowed_targets,
    const int min_poss_borderness,
    std::pair<int, int>& current_best_cell, int& current_best_borderness) {

    //cout << "call " << (++num_calls) << endl;
    //num_calls = 0;

    if (tile.empty()) {
        return 0;
    }

    const long long cost_base = 32;
    const long long cost_iter = 1;
    long long time_used = cost_base;

    auto iter = tile.begin();
    while (iter != tile.end()) {// && (*iter).borderness < current_best_borderness) {
        auto entry = (*iter);
        std::pair<int, int> probe = entry.pos_pair();
        iter++;
        time_used += cost_iter;
        //if (current_best_borderness != INT_MAX && rand() % 1024 < 256)
        //    continue;
        if (board->query_conflict_free(probe)) {
            allowed_targets.push_back(probe);
            /*current_best_cell = probe;
            current_best_borderness = entry.borderness;
            if (current_best_borderness <= min_poss_borderness) {
                break;
            }*/
        }
    }

    //if(num_calls)
    //cout << num_calls << endl;

    return time_used;
}
