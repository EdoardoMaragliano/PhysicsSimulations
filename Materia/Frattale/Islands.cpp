#include "Islands.hpp"
#include <vector>
#include <iostream>

using namespace std;

Islands::Islands(){
    ROW=10;
    COL=10;
}
int Islands::PC(int num){           
  if (num == ROW){num = 0;}
  else if (num == -1){num = ROW-1;}
  return num;
}

// A function to check if a given
// cell (row, col) can be included in DFS
int Islands::isSafe(vector<vector<int>> M, int row, int col, vector<vector<bool>> visited)
{
    // row number is in range, column
    // number is in range and value is 1
    // and not yet visited
    return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL) && (M[row][col] && !visited[row][col]);

}
 
// A utility function to do DFS for a
// 2D boolean matrix. It only considers
// the 8 neighbours as adjacent vertices
void Islands::DFS(vector<vector<int>> M, int row, int col, vector<vector<bool>> visited)
{
    // These arrays are used to get
    // row and column numbers of 4
    // neighbours of a given cell
    static int rowNbr[] = { -1, 0, 0, 1};
    static int colNbr[] = {  0,-1, 1, 0};
 
    // Mark this cell as visited
    visited[PC(row)][PC(col)] = true;
 
    // Recur for all connected neighbours
    for (int k = 0; k < 4; ++k)
        if (isSafe(M, PC(row + rowNbr[k]), PC(col + colNbr[k]), visited))
            DFS(M, PC(row + rowNbr[k]), PC(col + colNbr[k]), visited);
}
 
// The main function that returns
// count of islands in a given boolean
// 2D matrix
int Islands::countIslands(vector<vector<int>> M){

    //Update private members to matrix dimensions (square)
    ROW=M.size();
    COL=M[0].size();
    cout << "ROW="<<ROW<<endl<<"COL="<<COL<<endl;

    // Make a bool array to mark visited cells.
    // Initially all cells are unvisited
    vector<vector<bool>> visited;

    for (int i = 0; i < ROW; ++i) { 
        vector<bool> v_appo;
        for (int j = 0; j < COL; ++j) { 
            v_appo.push_back(false); 
        } 
    visited.push_back(v_appo); 
    } 
    cout <<"visited=NULL"<<endl;
    // Initialize count as 0 and
    // travese through the all cells of
    // given matrix
    int count = 0;
    for (int i = 0; i < ROW; ++i)
        for (int j = 0; j <COL; ++j)
 
            // If a cell with value 1 is not  
            // visited yet, then new island found
            // Visit all cells in this island.
            if (M[i][j] && !visited[i][j]) {
                
                DFS(M, i, j, visited);
                cout <<"DFS"<<endl;
                // and increment island count
                ++count;
            }
 
    return count;
}
