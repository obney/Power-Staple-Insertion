#ifndef TYPES_H
#define TYPES_H

#include <unordered_set>

struct Chip {
    int left;
    int bottom;
    int right;
    int top;
    int num_rows;
    int row_height;
    int site_width;
    int num_cell_types;
    int num_cells;
};

struct CellType {
    int index;
    int width;
    int height;
    std::unordered_set<int> pin_sites;
};

struct Cell {
    int index;
    int type_index;
    int x;
    int y;
    int max_disp;
};

struct CellOut {
    int index;
    int x;
    int y;
    bool flip;
};

struct Staple {
    int x;
    int y;
};

#endif