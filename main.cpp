#include "testing-celllist.h"
using namespace std;

//namespace CellList{

int main(){
CellList::CellList list;
array<double,3> c1 = {0.2,0.2,0.2};
array<double,3> c2 = {0.0,0.0,0.0};
array<double,3> c3 = {0.95,0.95,0.95};
vector<array<double,3>> coords;
coords.push_back(c1);
coords.push_back(c2);
coords.push_back(c3); 
double box[3] = {1,1,1};
double distance = list.distance(coords[0], coords[1]);
//cout << distance << "\n";
array<int,3> boxdim = list.BoxtoCells(box, 0.4);
cout << "box dim" << boxdim[0] << ";" << boxdim[1] << ";" << boxdim[2] << "\n";
map<array<int,3>, set<int>> cellgen = list.Cellgen(box, coords, 3, 0.25);
array<int,3> cin = {0,0,0};
//set<array<int,3>> neigh = list.findCellneighbours(box, cin, 0.4);
map<int, set<int>> celllist = list.CellListgen(box, coords, cellgen, 3, 0.25);

for (auto it = celllist.begin(); it != celllist.end(); it++){
cout << it->first << ": " ;
auto set = it->second;
for( auto s = set.begin(); s != set.end(); s++){
    cout << *s << "\n " ;
}

};
return 0;
};

//}// namespace CellList