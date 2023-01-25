#include "testing-celllist.h"

using namespace std;

namespace CellList{

// calculating the distance between two coordinate points with pythagoras theorem
double CellList::distance(point coord1, point coord2){
return std::sqrt(std::pow(coord2[0]-coord1[0],2) + std::pow(coord2[1]- coord1[1], 2)+ std::pow(coord2[2]- coord1[2], 2));
};

// calculate nr of cells in each box dimensions, taking ceiling of side length/cutoff (rounding up)
std::array<int,3> CellList::boxtoCells(boxdim box, const double cutoff){
int nrcells[3] = {int(ceil(box[0]/cutoff)),int(ceil(box[1]/cutoff)),int(ceil(box[2]/cutoff))};
std::array<int,3> nrcells_array {nrcells[0], nrcells[1], nrcells[2]};
// total number of cells in box (should be positive)
//int totcells = abs(nrcells[0]*nrcells[1]*nrcells[2]);
return nrcells_array;
};

// find which cell particle is in: calculate how many "cutoffs" from origin particle is in each direction
cellindex CellList::findCell(point coord, double cutoff){
cellindex index = {int(floor(coord[0]/cutoff)),int(floor(coord[1]/cutoff)),int(floor(coord[2]/cutoff))};
return index; // returning index of cell point is located in
};

// creating a map with cell index as key and value a set of all point-indices located within corresponding cell
std::map<cellindex, set<pointindex>> CellList::genCell(boxdim box, std::vector<point> coords, int nr_particles, double cutoff){
std::array<int,3> nrcells = boxtoCells(box, cutoff); // calculating nr of cells in each dimension of box
// assign each particle to cell
std::multimap<cellindex, pointindex> coords_in; // cell index + point index
for(int i = 0; i < nr_particles; i++){
// find which cell particle is in: calculate how many "cutoffs" particle is in each direction
cellindex index = findCell(coords[i], cutoff);
coords_in.insert({index, i});
}
// find all points that are in same cell and put in a set to make a map with specific keys
std::map<cellindex, set<pointindex>> cellPIndex; // cell index + list of all points in that cell
for (int x =0; x<nrcells[0]; x++){ // loop over each dimension to go through all cells
    for(int y=0; y<nrcells[1]; y++){
        for(int z=0; z<nrcells[2]; z++){
            std::set<pointindex> pts_index;
            cellindex c_index = {x,y,z}; // defining cell index
            //cout << x << ", " << y << ", " << z << "\n";
            auto pts = coords_in.equal_range(c_index); // searching multimap for all values where key is c_index
            if(pts.first!=pts.second){ // cell not empty
                for(auto it = pts.first; it!=pts.second; it++){
                    //cout << it->second << "\n";
                    pts_index.insert(it->second); // adding all point indices to set
                }
                cellPIndex.insert({c_index, pts_index}); // adding c_index as key and its corresponding set with all points located in the cell  
            }
            if(pts.first == pts.second){// cell empty
                cellPIndex.insert({c_index, pts_index});
            }
        }
    }
}
return cellPIndex;
};


// finding all neighbour cells and add to one set
std::set<cellindex> CellList::findCellneighbours(boxdim box, const cellindex& c_index, const double cutoff){
// construct output set
std::set<cellindex> neigh_cells;
// construct sets to loop over for respective coordinate to find all 26 neighbours
// neighbour cells are +/1 in each dim 
std::set<int> x {c_index[0], (c_index[0]+1), (c_index[0]-1)}; 
std::set<int> y {c_index[1], (c_index[1]+1), (c_index[1]-1)};
std::set<int> z {c_index[2], (c_index[2]+1), (c_index[2]-1)};
// loop over indices and add neightbour cells to neigh_cells set

for (auto _x : x)
{
    for (auto _y : y)
    {
        for (auto _z : z)
        {
            cellindex new_index = {_x, _y, _z};
            neigh_cells.insert(new_index);
        }
    }
}
return neigh_cells;
};


// generating the Cell List
std::map< pointindex, set<pointindex>> CellList::genCellList(boxdim box, std::vector<point> coords, std::map<cellindex, set<pointindex>> indexpts, int nr_particles, double cutoff){
std::map< pointindex, set<pointindex>> celllist;
std::array<int,3> nrcells = boxtoCells(box, cutoff); // calculate nr of cells in each dim
// loop over all points/particles
for (int i = 0; i < nr_particles; i++){
    cout << i << " pts1 \n";
    cellindex cell = findCell(coords[i], cutoff); // find cell index where point i is
    std::set<cellindex> cell_neigh = findCellneighbours(box, cell, cutoff); // finding all neightbouring cells
    std::set<pointindex> all_pts;
    std::set<pointindex> ghost_pts;
    std::map< point,pointindex> pts_coords;
    // add all points located in the neighbouring cells to a set
    for( auto c = cell_neigh.begin(); c != cell_neigh.end(); c++){
        auto cindex = indexpts.find(*c);
        if(cindex!=indexpts.end()){ // the index was in map - c is one of the cells
           auto pts = cindex->second;
            if(pts.empty()==0){ // checking that set is not empty/cell is not empty
                all_pts.insert(pts.begin(), pts.end());
            }  
        }
        if(cindex == indexpts.end()){ // cell index not found due to perioditicity - c is one of the "ghost" cells
            //cout << "ghostcell" << "\n";
            cellindex ghostcell = *c;
            //cout << ghostcell[0] << ghostcell[1] << ghostcell[2] << "\n";
            cellindex ghostindex = {(nrcells[0] + (ghostcell[0]%nrcells[0]))%nrcells[0], (nrcells[1] + (ghostcell[1]%nrcells[1]))%nrcells[1], (nrcells[2] + (ghostcell[2]%nrcells[2]))%nrcells[2]}; // with modulo find which cells the perioditicity corresponds to
            //cout << ghostindex[0] << ghostindex[1] << ghostindex[2] << "\n";
            auto cindex = indexpts.find(ghostindex);
            auto pts = cindex->second;
            if( pts.empty()==0){ // checking that set is not empty/cell is not empty
                std::array<double,3> cell_diff = {double(ghostcell[0] - ghostindex[0]), double(ghostcell[1]-ghostindex[1]), double(ghostcell[2]-ghostindex[2])};
                for(auto p = pts.begin(); p!=pts.end();p++){
                    auto ghostcoord = coords[*p];
                    cout << ghostcoord[0] << ghostcoord[1] << ghostcoord[2] << "\n";
                    point new_coord = {ghostcoord[0]+(cutoff*cell_diff[0]), ghostcoord[1]+(cutoff*cell_diff[1]), ghostcoord[2]+(cutoff*cell_diff[2])};
                    cout << new_coord[0] << new_coord[1] << new_coord[2] << "\n";
                    pts_coords.insert({new_coord, *p});
                }
            }
        }
    }
    for( auto pts = all_pts.begin(); pts != all_pts.end(); pts++){ // add the "normal points" to neighbouring points coordinate map
        pts_coords.insert({coords[*pts], *pts});
    }
    std::set<pointindex> inter_pts;
    // go through points from neighbouring cells to see which ones are within cutoff
    for( auto pts2 = pts_coords.begin(); pts2 != pts_coords.end(); pts2++){
        cout << i << pts2->second << "point inter \n";
        if( i != pts2->second){
            double d = distance(coords[i], pts2->first);
            cout << "dis: " << d << "\n";
            if( d <= cutoff){
                inter_pts.insert(pts2->second);
            }
        }
    }
    celllist.insert({i, inter_pts}); // inserting key point index and value a set with point indices to all points within cutoff
}
return celllist; // returning cell list as a map
};

std::set<pointindex> CellList::genPointList(int pointin, boxdim box, std::vector<point> coords, std::map<cellindex, set<pointindex>> indexpts, int nr_particles, double cutoff){
std::set<pointindex> ptslist; // set to save point indices in
cellindex cell = findCell(coords[pointin], cutoff); // find cell index where input point is
std::set<cellindex> cell_neigh = findCellneighbours(box, cell, cutoff); // finding all neightbouring cells
std::set<pointindex> all_pts;
// add all points located in the neighbouring cells to a set
for( auto c = cell_neigh.begin(); c != cell_neigh.end(); c++){
    auto pts = indexpts.find(*c)->second;
    all_pts.insert(pts.begin(), pts.end());
}
// go through points from neighbouring cells to see which ones are within cutoff to the input point
for( auto pts2 = all_pts.begin(); pts2 != all_pts.end(); pts2++){
    if( pointin != *pts2){
        double d = distance(coords[pointin], coords[*pts2]);
        if( d <= cutoff){
            ptslist.insert(*pts2);
        }
    }
}
return ptslist; // resturn set of indices to points within cutoff from the input point
};

// after movement update which cells each point is located in 
map<array<int,3>, set<int>> CellList::updatePositions(vector<array<double,3>> old_coords, map<array<int,3>, set<int>> old_cellptsin, vector<array<double,3>> new_coords, int nr_particles, double cutoff){
for( int i =0 ; i< nr_particles; i++){ // loop over points
    array<int,3> old_cell = findCell(old_coords[i], cutoff);
    array<int, 3> new_cell = findCell(new_coords[i], cutoff);
    if (old_cell != new_cell){ // check if point moved enough to have changed cell
        auto erase = old_cellptsin.find(old_cell)->second;
        erase.erase(i); // erase point index from old cell location
        old_cellptsin[old_cell] = erase;
        auto add = old_cellptsin.find(new_cell)->second;
        add.insert(i); // add point index to new cell location
        old_cellptsin[new_cell] = add;
    }
}
return old_cellptsin; //now technically this is the new_cellptsin (cell index + points indices located in respective cell)
};

// after movement of one point update which cell point is located in 
std::map<cellindex, std::set<pointindex>> CellList::updateOnePosition(pointindex pointin, std::vector<point> old_coords, map<cellindex, set<pointindex>> old_cellptsin, std::vector<point> new_coords, int nr_particles, double cutoff){
cellindex old_cell = findCell(old_coords[pointin], cutoff); // find which cell point was located in    
cellindex new_cell = findCell(new_coords[pointin], cutoff); // find which cell point is now located in
if (old_cell != new_cell){ // check in point moved enough to have changed cell
    auto erase = old_cellptsin.find(old_cell)->second;
    erase.erase(pointin); // erase point index from old cell location
    old_cellptsin[old_cell] = erase;
    auto add = old_cellptsin.find(new_cell)->second;
    add.insert(pointin); // add point index to new cell location
    old_cellptsin[new_cell] = add;
}
return old_cellptsin; //now technically this is the new_cellptsin (cell index + points indices located in respective cell)  
};

}// namespace CellList