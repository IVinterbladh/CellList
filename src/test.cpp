#include "celllist.h"


//namespace CellList{

int main(){
CellList::CellList list;
/*
std::array<double,3> c1 = {0.2,0.2,0.2};
std::array<double,3> c2 = {0.0,0.0,0.0};
std::array<double,3> c3 = {0.95,0.95,0.95};
std::vector<std::array<double,3>> coords;
coords.push_back(c1);
coords.push_back(c2);
coords.push_back(c3); 
std::array<double,3> box = {1,1,1};
double distance = list.distance(coords[0], coords[1]);
//cout << distance << "\n";
//std::array<int,3> boxdim = list.BoxtoCells(box, 0.4);
//cout << "box dim" << boxdim[0] << ";" << boxdim[1] << ";" << boxdim[2] << "\n";
std::map<std::array<int,3>, std::set<int>> cellgen = list.genCell(box, coords, 3, 0.25);
std::array<int,3> cin = {0,0,0};
//set<array<int,3>> neigh = list.findCellneighbours(box, cin, 0.4);
std::map<int, std::set<int>> celllist = list.genCellList(box, coords, cellgen, 3, 0.25);

for (auto it = celllist.begin(); it != celllist.end(); it++){
std::cout << it->first << ": " ;
auto set = it->second;
for( auto s = set.begin(); s != set.end(); s++){
    std::cout << *s << "\n " ;
}

}; */
return 0;
};

//}// namespace CellList
