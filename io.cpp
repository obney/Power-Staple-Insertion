#include "io.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

extern Chip chip;
extern vector<CellType> cell_types;
extern vector<Cell> cells;
extern vector<vector<int>> R;
extern vector<vector<int>> ori_grid;
extern vector<vector<int>> grid;
extern vector<vector<CellOut>> cell_output;
extern vector<Staple> staples;

const int seed_num = 7777; //12345

void read_input(const string& input_path) {
    ifstream input_file(input_path);
    string line;

    getline(input_file, line);
    istringstream iss(line);
    iss >> chip.left >> chip.bottom >> chip.right >> chip.top;

    getline(input_file, line);
    iss.clear();
    iss.str(line);
    iss >> chip.num_rows >> chip.row_height >> chip.site_width;

    getline(input_file, line);
    chip.num_cell_types = stoi(line);

    getline(input_file, line);
    chip.num_cells = stoi(line);

    ori_grid.assign(chip.num_rows, vector<int>((chip.right - chip.left) / chip.site_width, 0));
    grid.assign(chip.num_rows, vector<int>((chip.right - chip.left) / chip.site_width, 0));
    R.assign(chip.num_rows, vector<int>());

    cell_types.resize(chip.num_cell_types);
    for (int i = 0; i < chip.num_cell_types; i++) {
        getline(input_file, line);
        istringstream iss(line);
        int index, width, height;
        iss >> index >> width >> height;
        unordered_set<int> pin_sites;
        int pin;
        while (iss >> pin) {
            pin_sites.insert(pin);
        }
        cell_types[index] = {index, width, height, pin_sites};
    }

    srand(seed_num);
    cells.resize(chip.num_cells);
    for (int i = 0; i < chip.num_cells; i++) {
        getline(input_file, line);
        istringstream iss(line);
        int idx, type_idx, x, y, max_disp;
        iss >> idx >> type_idx >> x >> y >> max_disp;
        int row = (y - chip.bottom) / chip.row_height;
        R[row].push_back(idx);
        int rand_num = rand();
        if((chip.num_cells > 200000 || (chip.right - chip.left) > 700000)) {
            if(chip.num_cells > 230000 && rand_num < RAND_MAX / 2)
                max_disp = max(0, max_disp - chip.site_width);
            if(chip.num_cells <= 230000 && rand_num < RAND_MAX / 5)
                max_disp = max(0, max_disp - chip.site_width);
        }
        cells[idx] = {idx, type_idx, x, y, max_disp};
        for(auto& pin : cell_types[type_idx].pin_sites)
            ori_grid[row][(x - chip.left) / chip.site_width + pin] = 1;
    }

    input_file.close();

    for(int i = 0; i < chip.num_rows; i++) {
        sort(R[i].begin(), R[i].end(), [&](int a, int b) {
            return cells[a].x < cells[b].x;
        });
    }
}

void write_output(const string& output_path) {
    ofstream output_file(output_path);

    for(int i = 0; i < cell_output.size(); i++) {
        for(auto& cell_out : cell_output[i]) {
            output_file << cell_out.index << " " << cell_out.x << " " << cell_out.y << " " << (cell_out.flip ? 1 : 0) << '\n';
        }
    }
    for(auto& staple : staples) {
        output_file << staple.x << " " << staple.y << '\n';
    }
    output_file.close();
}