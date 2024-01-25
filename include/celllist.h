#include <array>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <iostream>
#include <math.h> 


namespace CellList{
 
    using point = std::array<double,3>;
    using pointindex = int;
    using boxdim = std::array<double,3>;
    using cellindex = std::array<int,3>;

class CellList{

private:
boxdim box; // dimeansion of box in 3D, boxdim
double cutoff; // provided value of cutoff for the particles, givencutoff
std::array<int,3> boxcelldim; // nr of cells in each dimension of box
std::map<cellindex, std::set<pointindex>> cellindexMap; // map of cells with indices to points located in respective cell

///
/// @brief Calculating number of cells in box depending on cutoff
///        
/// @details void, but boxcelldim parameter is a filled array with number of cells in each dimension
///
void updateBoxtoCells();

public:
///
/// @brief Setting up the given parameters
///        
/// @details void, box dimensions and cutoff set for the private variables and boxcelldim and cellindexMap constructed
///   
/// @param[in] boxdim: provided dimensions of box; height, length, width
/// @param[in] givencutoff: provided value used as cutoff, beyond cutoff no interaction between particles
/// @param[in] coords: the coordinates of each point in x,y,z
/// @param[in] nr_points: int of number of particles/point
/// 
template<typename Pvector> // typename Pvector representing list/vector/etc of the point coordinates
void set_upCellList(boxdim boxdim, const double givencutoff, Pvector coords, int nr_points);

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
template< typename P> //typename P representing type of point ´
double distance(P coord1, P coord2);

///
/// @brief Finding which cell a point is located in
///        
/// @details 
///   
/// @param[in] coord: the coordinates of point; x,y,z
///
/// @return array with the cell index
///
template<typename P> //typename P representing type of point ´
cellindex findCell(P coord);

///
/// @brief Generate map to find which points are located in each cell
///        
/// @details map with cell index as key and value a set of point-indices located in respective cell
///
/// @param[in] coords: the coordinates of each point in x,y,z
/// @param[in] nr_points: number of particles/points
///
/// @return map with an array of the cell index and as value a set of all point-indices 
///
template<typename Pvector> // typename Pvector representing list/vector/etc of the point coordinates
void generateCell(Pvector coords, int nr_points);


///
/// @brief Finding all neighbours to a cell 
///        
/// @details Finding all neightbours to one specific cell including those due to periodicity 
///     
/// @param[in] c_index: the index of the cell whose neightbours are to be found
///
/// @return a set of all the indices of the neightbour cells
///
std::set<cellindex> findCellneighbours(const cellindex& c_index); // change to vector

///
/// @brief Compute List with all points located within the cutoff from one specific point
///        
/// @details 
/// 
/// @param[in] pointin: index of point for which interacting points should be found    
/// @param[in] coords: the coordinates of each point in x,y,z
/// @param[in] indexpts: a map providing cell as key and value which points are located within it
/// @param[in] nr_points: number of particles/points
///
/// @return a set containing all indices to points interacting with the input point
///
template<typename Pvector>
std::set<pointindex> generatePointList(const pointindex& pointin, Pvector coords, int nr_points);

///
/// @brief Computing the Cell List with all coordinate combinations located within the cutoff
///        
/// @details 
///     
/// @param[in] coords: the coordinates of each point in x,y,z
/// @param[in] indexpts: a map providing cell as key and value which points are located within it
/// @param[in] nr_points: number of particles/points
///
/// @return a Cell List as a map, key is a point index and the value is a set containing all indices to points interacting with the key point
///
template<typename Pvector>
std::map< pointindex, std::set<pointindex>> generateCellList(Pvector coords, int nr_points);

///
/// @brief After movement of one point- update position in cell
///        
/// @details 
///         
/// @param[in] pointin: index of point which has changed position 
/// @param[in] old_pcoord: coordinate values for old location of point
/// @param[in] new_pcoord: coordinate values for new location of point
///
/// @return updated map with cells as key and their respective values being a set of points located in them -> new_cellptsin
///
template<typename P>
void updateOnePosition(const pointindex& pointin, P old_pcoord, P new_pcoord);

///
/// @brief After movement (possibly all points changed)- update point positions in respective cell
///        
/// @details 
///          
/// @param[in] old_coords: the coordinates of each point in x,y,z berfore movement
/// @param[in] new_coords: the coordinates of each point in x,y,z after movement
/// @param[in] nr_points: number of particles/points
///
/// @return updated map with cells as key and their respective values being a set of points located in them -> new_cellptsin
///
template<typename Pvector>
void updatePositions(Pvector old_coords, Pvector new_coords, int nr_points);


};

}// namespace CellList