#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <array>
#include <functional>
#include <cmath>
#include <limits.h>
#include <omp.h>
#include "types.h"
#include "io.h"

using namespace std;

struct ArrayHasher {
    std::size_t operator()(const array<int, 7>& arr) const {
        std::size_t hash = 0;
        for (int i : arr) {
            hash ^= std::hash<int>{}(i) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

int thread_num = 8;
const int DP_row_size = 3;
Chip chip;
vector<CellType> cell_types;
vector<Cell> cells;
vector<vector<int>> R;
vector<vector<int>> ori_grid;
vector<vector<int>> grid;
vector<vector<CellOut>> cell_output(8, vector<CellOut>());
vector<Staple> staples;
vector<Staple> type1_staples;
vector<Staple> type2_staples;

int estimate(array<int, 7>& arg, int top_row, bool is_first, int len) {
    int added = 0, cur = arg[0] - 1;
    
    bool space00 = false, space0 = false, space1 = false, space2 = false, space3 = false;
    if(!is_first) {
        space00 = grid[top_row + 2][cur] == 0;
        space0 = grid[top_row + 1][cur] == 0;
    }
    if(len >= 1) space1 = (arg[1] == -1 || cell_types[cells[R[top_row][arg[1]]].type_index].pin_sites.find(arg[2] - 1) == cell_types[cells[R[top_row][arg[1]]].type_index].pin_sites.end());
    if(len >= 2) space2 = (arg[3] == -1 || cell_types[cells[R[top_row - 1][arg[3]]].type_index].pin_sites.find(arg[4] - 1) == cell_types[cells[R[top_row - 1][arg[3]]].type_index].pin_sites.end());
    if(len >= 3) space3 = (arg[5] == -1 || cell_types[cells[R[top_row - 2][arg[5]]].type_index].pin_sites.find(arg[6] - 1) == cell_types[cells[R[top_row - 2][arg[5]]].type_index].pin_sites.end());

    if(space0 && space1) added+=1;
    if(space1 && space2) added+=1;
    if(space2 && space3) added+=1;
    if(!space00 && space0 && !space1) added -= 1;
    if(!space0 && space1 && !space2) added -=1;
    if(!space1 && space2 && !space3) added -=1;
    
    
    return added;
}

//idx      0   1   2   3   4   5   6     
//arg = {cur, s0, l0, s1, l1, s2, l2};
int dfs(array<int, 7>& arg, int top_row, unordered_map<array<int, 7>, int, ArrayHasher>& memo, unordered_map<array<int, 7>, array<int, 7>, ArrayHasher>& next, bool is_first, int len) {
    if(memo.find(arg) != memo.end()){
        return memo[arg];
    }

    // cout << arg[0] << " " << arg[1] << " " << arg[2] << " " << arg[3] << " " << arg[4] << " "<< arg[5] << " " << arg[6] << '\n';
    if(chip.left + arg[0] * chip.site_width >= chip.right){
        if(arg[1] == ((int)R[top_row].size() - 1) && (len <= 1 || arg[3] == ((int)R[top_row - 1].size() - 1)) && (len <= 2 || arg[5] == ((int)R[top_row - 2].size() - 1))) {
            for(int i = 0; i < len; i++)
                if(arg[i * 2 + 1] != -1 && arg[i * 2 + 2] * chip.site_width < cell_types[cells[R[top_row - i][arg[i * 2 + 1]]].type_index].width) return memo[arg] = INT_MIN;
            return memo[arg] = 0;
        }
        return memo[arg] = INT_MIN;
    }
    

    int benefit = INT_MIN;
    array<array<bool, 3>, 2> choice = {{{1, 1, 1}, {0, 0, 0}}};

    for(int i = 0; i < len; i++) {
        if(arg[i * 2 + 1] == (int)R[top_row - i].size() - 1) {
            choice[0][i] = true;
            choice[1][i] = false;
            continue;
        }
        Cell& cur_cell = cells[R[top_row - i][arg[i * 2 + 1]]];
        Cell& next_cell = cells[R[top_row - i][arg[i * 2 + 1] + 1]];
        //defer
        choice[0][i] = ((next_cell.x + next_cell.max_disp) > chip.left + arg[0] * chip.site_width) ? true : false;
        //put
        choice[1][i] = (((next_cell.x - next_cell.max_disp) <= chip.left + arg[0] * chip.site_width) && ((next_cell.x + next_cell.max_disp) >= chip.left + arg[0] * chip.site_width) && (arg[i * 2 + 1] == -1 || cell_types[cur_cell.type_index].width <= arg[i * 2 + 2] * chip.site_width)) ? true : false;
    }


    for(int i = 0; i < 8; i++) {
        // cout << ((i & 1) > 0) << " " << ((i & 2) > 0) << " " << ((i & 4) > 0) << '\n';
        // cout << choice[(i & 1) > 0][0] << " " << choice[(i & 2) > 0][1] << " " << choice[(i & 4) > 0][2] << '\n';
        if(choice[(i & 1) > 0][0] && choice[(i & 2) > 0][1] && choice[(i & 4) > 0][2]) {
            array<int, 7> new_arg = arg;
            new_arg[0] += 1;
            bool out_of_bound = false;
            for(int j = 0; j < len; j++) {
                new_arg[j * 2 + 1] = (i & (1 << j)) ? arg[j * 2 + 1] + 1 : arg[j * 2 + 1];
                new_arg[j * 2 + 2] = (i & (1 << j)) ? 1 : arg[j * 2 + 2] + 1;
                if(new_arg[j * 2 + 1] >= (int)R[top_row - j].size()) {
                    out_of_bound = true;
                    break;
                }
            }
            
            if(out_of_bound) continue;
            
            int re = dfs(new_arg, top_row, memo, next, is_first, len);
            if(re == INT_MIN) continue;
            int added = estimate(new_arg, top_row, is_first, len);
            if(re + added > benefit) {
                benefit = re + added;
                next[arg] = new_arg;
            }
        }
    }
    return memo[arg] = benefit;
}

void update_cell_position(unordered_map<array<int, 7>, array<int, 7>, ArrayHasher>& next, int top_row, int thread_id) {
    array<int, 7> key = {0, -1, 0, -1, 0, -1, 0};
    while(1) {
        //cout << key[0] << " " << key[1] << " " << key[2] << " " << key[3] << " " << key[4] << " " << key[5] << " " << key[6] << '\n';
        for(int i = 0; i < 3; i++) {
            if(key[i * 2 + 1] >= 0 && key[i * 2 + 2] == 1) {
                for(auto& pin: cell_types[cells[R[top_row - i][key[i * 2 + 1]]].type_index].pin_sites) {
                    grid[top_row - i][key[0] - 1 + pin] = 1;
                }
                cell_output[thread_id].push_back(CellOut{R[top_row - i][key[i * 2 + 1]], chip.left + (key[0] - 1) * chip.site_width, chip.bottom + (top_row - i) * chip.row_height, false});
            }
        }

        if(next.find(key) == next.end()) break;
        key = next[key];
    }
}


bool is_valid(int row, int col) {
    if(row - 1 < 0 || row > grid.size() || col < 0 || col > grid[0].size() || grid[row][col] != 0 || grid[row - 1][col] != 0) return false;
    if(((row + 1) < grid.size()) && ((col - 1) >= 0) && (grid[row + 1][col] <= 1) && (grid[row][col - 1] <= 1) && (grid[row + 1][col - 1] == 3))
        return false;
    grid[row - 1][col] = 3;
    if(((row - 2) >= 0) && ((col - 1) >= 0) && (grid[row - 1][col - 1] <= 1) && (grid[row - 2][col - 1] == 2)) {
        if(!is_valid(row - 2, col)) {
            grid[row - 1][col] = 0;
            return false;
        }
    }
    grid[row - 1][col] = 0;
    return true;
}

void greedy_put_staple() {
    for(int col = 0; col < grid[0].size(); col++) {
        int row = grid.size() - 1;
        while(row > 0) {
            if(is_valid(row, col)) {
                grid[row][col] = 2;
                grid[row - 1][col] = 3;
                if(((grid.size() - row) & 1) == 1) {
                    type1_staples.push_back(Staple{chip.left + col * chip.site_width, chip.bottom + (row - 1) * chip.row_height});
                }
                else {
                    type2_staples.push_back(Staple{chip.left + col * chip.site_width, chip.bottom + (row - 1) * chip.row_height});
                }
                row -= 2;
            }
            else row--;
        }
    }
}

int main(int argc, char* argv[]) {
    // auto start_time = chrono::high_resolution_clock::now();

    string input_path = argv[1];
    string output_path = argv[2];

    read_input(input_path);

    if(chip.num_rows <= 80 && (chip.right - chip.left) <= 120000) thread_num = 1;

    int block_size = ceil((double)chip.num_rows / thread_num);
    block_size += ((-block_size + 99999999 + (DP_row_size == 3 ? 0 : 1)) % (DP_row_size == 3 ? 3 : 2));

    //printGrid(ori_grid);

    #pragma omp parallel for num_threads(thread_num)
    for(int i = 0; i < thread_num; i++) {
        // cout << chip.num_rows - 1 - i * block_size << '\n';
        bool is_first = true;
        for(int j = chip.num_rows - 1 - i * block_size; j >= max(0, chip.num_rows - (i + 1) * block_size); j -= (DP_row_size == 3 ? 3 : 2)) {
            unordered_map<array<int, 7>, int, ArrayHasher> memo;
            unordered_map<array<int, 7>, array<int, 7>, ArrayHasher> next;
            array<int, 7> arg = {0, -1, 0, -1, 0, -1, 0};
            int len = 3;
            if(DP_row_size == 2) len = 2;
            if(j == 1) len = 2;
            else if(j == 0) len = 1;
            int benefit = dfs(arg, j, memo, next, is_first, len);
            // cout<< "Benefit: " << benefit << '\n';
            update_cell_position(next, j, i);
            is_first = false;
        }
    }
    
    greedy_put_staple();
    // cout << "Staple size: " <<  type1_staples.size() +  type2_staples.size() << '\n';
    // cout << "Type 1 staple size: " << type1_staples.size() << '\n';
    // cout << "Type 2 staple size: " << type2_staples.size() << '\n';
    if(type1_staples.size() > type2_staples.size()) swap(type1_staples, type2_staples);
    int maximun_staple_size = floor((double)type1_staples.size() * 1.1);
    if(type2_staples.size() > maximun_staple_size) {
        // cout << "imbalaced staple size" << '\n';
        int delete_size = type2_staples.size() - maximun_staple_size;
        while(delete_size > 0) {
            type2_staples.pop_back();
            delete_size--;
        }
        // cout << "Updated Staple size: " <<  type1_staples.size() +  type2_staples.size() << '\n';
    }
    // cout << "Ratio: " << type2_staples.size() / (double)type1_staples.size() << '\n';

    for(int i = 0; i < type1_staples.size(); i++) {
        staples.push_back(type1_staples[i]);
    }
    for(int i = 0; i < type2_staples.size(); i++) {
        staples.push_back(type2_staples[i]);
    }

    write_output(output_path);

    // auto end_time = chrono::high_resolution_clock::now();
    // auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    // cout << "Time taken: " << duration.count() << " ms" << '\n';

    //printGrid(grid);
    return 0;
}

