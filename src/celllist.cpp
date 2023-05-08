#include "celllist.h"


namespace CellList{
// calculate nr of cells in each box dimensions, taking ceiling of side length/cutoff (rounding up)
void CellList::updateBoxtoCells(){
boxcelldim = {int(ceil(box[0]/cutoff)),int(ceil(box[1]/cutoff)),int(ceil(box[2]/cutoff))};
};

// setting up values to private variables and constructing box dimensions
template<typename Pvector>
void CellList::set_upCellList(boxdim boxdim, const double givencutoff, Pvector coords, int nr_points){
    cutoff = givencutoff;
    box = boxdim;
    updateBoxtoCells();
    generateCell(coords, nr_points);
}

// calculating the distance between two coordinate points with pythagoras theorem
template< typename P>
double CellList::distance(P coord1, P coord2){
return std::sqrt(std::pow(coord2[0]-coord1[0],2) + std::pow(coord2[1]- coord1[1], 2)+ std::pow(coord2[2]- coord1[2], 2));
};

// find which cell particle is in: calculate how many "cutoffs" from origin particle is in each direction
template<typename P>
cellindex CellList::findCell(P coord){
cellindex index = {int(floor(coord[0]/cutoff)),int(floor(coord[1]/cutoff)),int(floor(coord[2]/cutoff))};
return index; // returning index of cell where point is located
};

// creating a map with cell index as key and value a set of all point-indices located within corresponding cell
template<typename Pvector>
void CellList::generateCell(Pvector coords, int nr_points){
// assign each particle to cell
std::multimap<cellindex, pointindex> coords_in; // cell index + point index
for(int i = 0; i < nr_points; i++){
// find which cell particle is in: calculate how many "cutoffs" particle is in each direction
cellindex index = findCell(coords[i]);
coords_in.insert({index, i});
}
// find all points that are in same cell and put in a set to make a map with specific keys
//std::map<cellindex, set<pointindex>> cellPIndex; // cell index + list of all points in that cell
for (int x =0; x<boxcelldim[0]; x++){ // loop over each dimension to go through all cells
    for(int y=0; y<boxcelldim[1]; y++){
        for(int z=0; z<boxcelldim[2]; z++){
            std::set<pointindex> pts_index;
            cellindex c_index = {x,y,z}; // defining cell index
            //cout << x << ", " << y << ", " << z << "\n";
            auto pts = coords_in.equal_range(c_index); // searching multimap for all values where key is c_index
            if(pts.first!=pts.second){ // cell not empty
                for(auto it = pts.first; it!=pts.second; it++){
                    //cout << it->second << "\n";
                    pts_index.insert(it->second); // adding all point indices to set
                }
                cellindexMap.insert({c_index, pts_index}); // adding c_index as key and its corresponding set with all points located in the cell  
            }
            if(pts.first == pts.second){// cell empty
                cellindexMap.insert({c_index, pts_index});
            }
        }
    }
}
};


// finding all neighbour cells and add to one set
std::set<cellindex> CellList::findCellneighbours(const cellindex& c_index){
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

template<typename Pvector>
std::set<pointindex> CellList::generatePointList(const pointindex& pointin, Pvector coords, int nr_points){
std::set<pointindex> ptslist; // set to save point indices in
cellindex cell = findCell(coords[pointin]); // find cell index where input point is
std::set<cellindex> cell_neigh = findCellneighbours(cell); // finding all neightbouring cells
std::set<pointindex> all_pts;
std::map< point,pointindex> pts_coords;
// add all points located in the neighbouring cells to a set
for( auto c = cell_neigh.begin(); c != cell_neigh.end(); c++){
    auto cindex = cellindexMap.find(*c);
    if(cindex!=cellindexMap.end()){ // the index was in map - c is one of the cells
        auto pts = cindex->second;
        if(pts.empty()==0){ // checking that set is not empty/cell is not empty
            all_pts.insert(pts.begin(), pts.end());
        }
    }
    if(cindex == cellindexMap.end()){ // cell index not found due to perioditicity - c is one of the "ghost" cells
    cellindex ghostcell = *c;
    cellindex ghostindex = {(boxcelldim[0] + (ghostcell[0]%boxcelldim[0]))%boxcelldim[0], (boxcelldim[1] + (ghostcell[1]%boxcelldim[1]))%boxcelldim[1], (boxcelldim[2] + (ghostcell[2]%boxcelldim[2]))%boxcelldim[2]}; // with modulo find which cells the perioditicity corresponds to
    auto cindex = cellindexMap.find(ghostindex);
    auto pts = cindex->second;
    if( pts.empty()==0){ // checking that set is not empty/cell is not empty
        std::array<double,3> cell_diff = {double(ghostcell[0] - ghostindex[0]), double(ghostcell[1]-ghostindex[1]), double(ghostcell[2]-ghostindex[2])};
        for(auto p = pts.begin(); p!=pts.end();p++){
        auto ghostcoord = coords[*p];
        point new_coord = {ghostcoord[0]+(cutoff*cell_diff[0]), ghostcoord[1]+(cutoff*cell_diff[1]), ghostcoord[2]+(cutoff*cell_diff[2])};
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
        if( pointin != pts2->second){
            double d = distance(coords[pointin], pts2->first);
            cout << "dis: " << d << "\n";
            if( d <= cutoff){
                inter_pts.insert(pts2->second);
            }
        }
    }
return ptslist; // resturn set of indices to points within cutoff from the input point
};

// generating the Cell List
template<typename Pvector>
std::map< pointindex, std::set<pointindex>> CellList::generateCellList(Pvector coords, int nr_points){
std::map< pointindex, std::set<pointindex>> celllist;
// loop over all points/particles
for (int i = 0; i < nr_points; i++){
    std::set<pointindex> pointlist = generatePointList(i, coords, nr_points);
    celllist.insert({i, pointlist}); // inserting key point index and value a set with point indices to all points within cutoff
}
return celllist; // returning cell list as a map
};

// after movement of one point update which cell point is located in 
template<typename P>
void CellList::updateOnePosition(const pointindex& pointin, P old_pcoord, P new_pcoord){
cellindex old_cell = findCell(old_pcoord); // find which cell point was located in    
cellindex new_cell = findCell(new_pcoord); // find which cell point is now located in
if (old_cell != new_cell){ // check if point moved enough to have changed cell
    auto erase = cellindexMap.find(old_cell)->second;
    erase.erase(pointin); // erase point index from old cell location
    cellindexMap[old_cell] = erase;
    auto add = cellindexMap.find(new_cell)->second;
    add.insert(pointin); // add point index to new cell location
    cellindexMap[new_cell] = add;
}//now technically this is the new_cellptsin (cell index + points indices located in respective cell)  
};

// after movement update which cells each point is located in
template<typename Pvector>
void CellList::updatePositions(Pvector old_coords, Pvector new_coords, int nr_points){
for( int i =0 ; i< nr_points; i++){ // loop over points
    updateOnePosition(i, old_coords[i], new_coords[i], nr_points); // calling function updating one position for every particle
}
};
}// namespace CellList
