#include <array>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iostream>
#include <math.h> 

using namespace std;

namespace CellList{

    using vec3 = array<double,3>;
    using CellCoord = array<int,3>;

class CellList{
private:

public:
//int nr_particles = 10;

///
/// @brief Calculating the distance 
///        
/// @details Calculate distance between two coordinates/points returning a double
///   
/// @param[in] coord1: double containing coordinates for a point
/// @param[in] coord2: double containing coordinates for a point
///
/// @return the distance as a double
///
double distance(array<double,3> coord1, array<double,3> coord2);

///
/// @brief Calculating number of cells in box depending on cutoff
///        
/// @details 
///   
/// @param[in] boxdim: dimensions of box; height, length, width
/// @param[in] cutoff: value used as cutoff, beyond cutoff no interaction between particles
///
/// @return array with number of cells in each dimension
///
array<int,3> BoxtoCells(double boxdim[3], double cutoff);

///
/// @brief Finding which cell a point is located in
///        
/// @details 
///   
/// @param[in] coord: the coordinates of point; x,y,z
/// @param[in] cutoff: value used as cutoff, beyond cutoff no interaction between particles
///
/// @return array with the cell index
///
array<int,3> findCell(array<double,3> coord, double cutoff);

///
/// @brief Generate map to find which points are located in each cell
///        
/// @details map with cell index as key and value a set of point-indices located in respective cell
///
/// @param[in] boxdim: the dimension of the original box       
/// @param[in] coords: vector of doubles containing the coordinates for respective point
/// @param[in] nr_particles: number of particles/points
/// @param[in] cutoff: the value of the cutoff
///
/// @return map with an array of the cell index and as value a set of all point-indices 
///
map<array<int,3>, set<int>> Cellgen(double boxdim[3], vector<array<double,3>> coords, int nr_particles, double cutoff);


///
/// @brief Finding all neighbours to a cell 
///        
/// @details Finding all neightbours to one specific cell including those due to periodicity 
///     
/// @param[in] boxdim: the dimension of the original box 
/// @param[in] c_index: the index of the cell whose neightbours are to be found
///
/// @return a set of all the indices of the neightbour cells
///
set<array<int,3>> findCellneighbours(double boxdim[3], const array<int,3>& c_index, double cutoff); // change to vector

///
/// @brief Computing the Cell List with all coordinate combinations located within the cutoff
///        
/// @details 
///     
/// @param[in] boxdim: the dimension of the original box 
/// @param[in] coords: array of doubles containing the coordinates for respective point and the cell number it is located in
/// @param[in] indexpts: a map providing cell as key and value which points are located within it
/// @param[in] nr_particles: number of particles/points
/// @param[in] cutoff: the value of the cutoff
///
/// @return a Cell List as a map, key is a point index and the value is a set containing all indices to points interacting with the key point
///
map< int, set<int>> CellListgen(double boxdim[3], vector<array<double,3>> coords, map<array<int,3>, set<int>> indexpts, int nr_particles, double cutoff);

///
/// @brief Compute List with all points located within the cutoff from one specific point
///        
/// @details 
///     
/// @param[in] boxdim: the dimension of the original box 
/// @param[in] coords: array of doubles containing the coordinates for respective point and the cell number it is located in
/// @param[in] indexpts: a map providing cell as key and value which points are located within it
/// @param[in] nr_particles: number of particles/points
/// @param[in] cutoff: the value of the cutoff
///
/// @return a set containing all indices to points interacting with the input point
///
set<int> pointListgen(int point, double boxdim[3], vector<array<double,3>> coords, map<array<int,3>, set<int>> indexpts, int nr_particles, double cutoff);


///
/// @brief After movement (possibly all points changed)- update point positions in respective cell
///        
/// @details 
///          
/// @param[in] old_coords: array of doubles containing the coordinates for respective point before movement
/// @param[in] old_cellptsin: map of before movement cell containing points indices 
/// @param[in] new_coords: array of doubles containing the coordinates for respective point after movement
/// @param[in] nr_particles: number of particles/points
/// @param[in] cutoff: the value of the cutoff
///
/// @return updated map with cells as key and their respective values being a set of points located in them -> new_cellptsin
///
map<array<int,3>, set<int>> updatePositions(vector<array<double,3>> old_coords, map<array<int,3>, set<int>> old_cellptsin, vector<array<double,3>> new_coords, int nr_particles, double cutoff);

///
/// @brief After movement of one point- update position in cell
///        
/// @details 
///         
/// @param[in] point: index of point which has changed position 
/// @param[in] old_coords: array of doubles containing the coordinates for respective point before movement
/// @param[in] old_cellptsin: map of before movement cell containing points indices 
/// @param[in] new_coords: array of doubles containing the coordinates for respective point after movement
/// @param[in] nr_particles: number of particles/points
/// @param[in] cutoff: the value of the cutoff
///
/// @return updated map with cells as key and their respective values being a set of points located in them -> new_cellptsin
///
map<array<int,3>, set<int>> updateOnePosition(int point, vector<array<double,3>> old_coords, map<array<int,3>, set<int>> old_cellptsin, vector<array<double,3>> new_coords, int nr_particles, double cutoff);
};

}// namespace CellList