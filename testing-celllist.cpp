#include "testing-celllist.h"

using namespace std;

namespace CellList{

// calculating the distance between two coordinate points with pythagoras theorem
double CellList::distance(array<double,3> coord1, array<double,3> coord2){
return sqrt(pow(coord2[0]-coord1[0],2) + pow(coord2[1]- coord1[1], 2)+ pow(coord2[2]- coord1[2], 2));
};

// calculate nr of cells in each box dimensions, taking ceiling of side length/cutoff (rounding up)
array<int,3> CellList::BoxtoCells(double boxdim[3], double cutoff){
int nrcells[3] = {int(ceil(boxdim[0]/cutoff)),int(ceil(boxdim[1]/cutoff)),int(ceil(boxdim[2]/cutoff))};
array<int,3> nrcells_array {nrcells[0], nrcells[1], nrcells[2]};
// total number of cells in box (should be positive)
//int totcells = abs(nrcells[0]*nrcells[1]*nrcells[2]);
return nrcells_array;
};

// find which cell particle is in: calculate how many "cutoffs" from origin particle is in each direction
array<int,3> CellList::findCell(array<double,3> coord, double cutoff){
array<int,3> index = {int(floor(coord[0]/cutoff)),int(floor(coord[1]/cutoff)),int(floor(coord[2]/cutoff))};
return index; // returning index of cell point is located in
};

// creating a map with cell index as key and value a set of all point-indices located within corresponding cell
map<array<int,3>, set<int>> CellList::Cellgen(double boxdim[3], vector<array<double,3>> coords, int nr_particles, double cutoff){
array<int,3> nrcells = BoxtoCells(boxdim, cutoff); // calculating nr of cells in each dimension of box
// assign each particle to cell
multimap<array<int,3>, int> coords_in; // cell index + point index
for(int i = 0; i < nr_particles; i++){
// find which cell particle is in: calculate how many "cutoffs" particle is in each direction
array<int,3> index = findCell(coords[i], cutoff);
coords_in.insert({index, i});
}
// find all points that are in same cell and put in a set to make a map with specific keys
map<array<int,3>, set<int>> cellindex; // cell index + list of all points in that cell
for (int x =0; x<nrcells[0]; x++){ // loop over each dimension to go through all cells
    for(int y=0; y<nrcells[1]; y++){
        for(int z=0; z<nrcells[2]; z++){
            set<int> pts_index;
            array<int,3> c_index = {x,y,z}; // defining cell index
            //cout << x << ", " << y << ", " << z << "\n";
            auto pts = coords_in.equal_range(c_index); // searching multimap for all values where key is c_index
            if(pts.first!=pts.second){ // cell not empty
                for(auto it = pts.first; it!=pts.second; it++){
                    //cout << it->second << "\n";
                    pts_index.insert(it->second); // adding all point indices to set
                }
                cellindex.insert({c_index, pts_index}); // adding c_index as key and its corresponding set with all points located in the cell  
            }
            if(pts.first == pts.second){// cell empty
                cellindex.insert({c_index, pts_index});
            }
        }
    }
}
return cellindex;
};


// finding all neighbour cells and add to one set
set<array<int,3>> CellList::findCellneighbours(double boxdim[3], const array<int,3>& c_index, const double cutoff){
// construct output set
set<array<int,3>> neigh_cells;
// construct sets to loop over for respective coordinate to find all 26 neighbours
// neighbour cells are +/1 in each dim 
set<int> x {c_index[0], (c_index[0]+1), (c_index[0]-1)}; 
set<int> y {c_index[1], (c_index[1]+1), (c_index[1]-1)};
set<int> z {c_index[2], (c_index[2]+1), (c_index[2]-1)};
// loop over indices and add neightbour cells to neigh_cells set

for (auto _x : x)
{
    for (auto _y : y)
    {
        for (auto _z : z)
        {
            array<int, 3> new_index = {_x, _y, _z};
            neigh_cells.insert(new_index);
        }
    }
}

for(auto xit = x.begin(); xit!=x.end(); ++xit){
    for(auto yit = y.begin(); yit!=y.end(); ++yit){
        for(auto zit = z.begin(); zit!=z.end(); ++zit){
            //cout << "c index" << *xit << *yit << *zit << "\n";
            array<int,3> new_index = {*xit,*yit,*zit};
            neigh_cells.insert(new_index);
        }
    }
}
return neigh_cells;
};



// generating the Cell List
map< int, set<int>> CellList::CellListgen(double boxdim[3], vector<array<double,3>> coords, map<array<int,3>, set<int>> indexpts, int nr_particles, double cutoff){
map< int, set<int>> celllist;
array<int,3> nrcells = BoxtoCells(boxdim, cutoff); // calculate nr of cells in each dim
// loop over all points/particles
for (int i = 0; i < nr_particles; i++){
    cout << i << " pts1 \n";
    array<int,3> cell = findCell(coords[i], cutoff); // find cell index where point i is
    set<array<int,3>> cell_neigh = findCellneighbours(boxdim, cell, cutoff); // finding all neightbouring cells
    set<int> all_pts;
    set<int> ghost_pts;
    map< array<double,3>,int> pts_coords;
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
            array<int,3> ghostcell = *c;
            //cout << ghostcell[0] << ghostcell[1] << ghostcell[2] << "\n";
            array<int,3> ghostindex = {(nrcells[0] + (ghostcell[0]%nrcells[0]))%nrcells[0], (nrcells[1] + (ghostcell[1]%nrcells[1]))%nrcells[1], (nrcells[2] + (ghostcell[2]%nrcells[2]))%nrcells[2]}; // with modulo find which cells the perioditicity corresponds to
            //cout << ghostindex[0] << ghostindex[1] << ghostindex[2] << "\n";
            auto cindex = indexpts.find(ghostindex);
            auto pts = cindex->second;
            if( pts.empty()==0){ // checking that set is not empty/cell is not empty
                array<double,3> cell_diff = {double(ghostcell[0] - ghostindex[0]), double(ghostcell[1]-ghostindex[1]), double(ghostcell[2]-ghostindex[2])};
                for(auto p = pts.begin(); p!=pts.end();p++){
                    auto ghostcoord = coords[*p];
                    cout << ghostcoord[0] << ghostcoord[1] << ghostcoord[2] << "\n";
                    array<double,3> new_coord = {ghostcoord[0]+(cutoff*cell_diff[0]), ghostcoord[1]+(cutoff*cell_diff[1]), ghostcoord[2]+(cutoff*cell_diff[2])};
                    cout << new_coord[0] << new_coord[1] << new_coord[2] << "\n";
                    pts_coords.insert({new_coord, *p});
                }
            }
        }
    }
    for( auto pts = all_pts.begin(); pts != all_pts.end(); pts++){ // add the "normal points" to neighbouring points coordinate map
        pts_coords.insert({coords[*pts], *pts});
    }
    set<int> inter_pts;
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

set<int> CellList::pointListgen(int point, double boxdim[3], vector<array<double,3>> coords, map<array<int,3>, set<int>> indexpts, int nr_particles, double cutoff){
set<int> ptslist; // set to save point indices in
array<int,3> cell = findCell(coords[point], cutoff); // find cell index where input point is
set<array<int,3>> cell_neigh = findCellneighbours(boxdim, cell, cutoff); // finding all neightbouring cells
set<int> all_pts;
// add all points located in the neighbouring cells to a set
for( auto c = cell_neigh.begin(); c != cell_neigh.end(); c++){
    auto pts = indexpts.find(*c)->second;
    all_pts.insert(pts.begin(), pts.end());
}
// go through points from neighbouring cells to see which ones are within cutoff to the input point
for( auto pts2 = all_pts.begin(); pts2 != all_pts.end(); pts2++){
    if( point != *pts2){
        double d = distance(coords[point], coords[*pts2]);
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
map<array<int,3>, set<int>> CellList::updateOnePosition(int point, vector<array<double,3>> old_coords, map<array<int,3>, set<int>> old_cellptsin, vector<array<double,3>> new_coords, int nr_particles, double cutoff){
array<int,3> old_cell = findCell(old_coords[point], cutoff); // find which cell point was located in    
array<int, 3> new_cell = findCell(new_coords[point], cutoff); // find which cell point is now located in
if (old_cell != new_cell){ // check in point moved enough to have changed cell
    auto erase = old_cellptsin.find(old_cell)->second;
    erase.erase(point); // erase point index from old cell location
    old_cellptsin[old_cell] = erase;
    auto add = old_cellptsin.find(new_cell)->second;
    add.insert(point); // add point index to new cell location
    old_cellptsin[new_cell] = add;
}
return old_cellptsin; //now technically this is the new_cellptsin (cell index + points indices located in respective cell)  
};

}// namespace CellList