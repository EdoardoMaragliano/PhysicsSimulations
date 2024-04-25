#ifndef _ISLANDS_
#define _ISLANDS_

#include <vector>
#include <iostream>


using namespace std;

class Islands{
public:
    Islands();
    // The main function that returns
    // count of islands in a given boolean
    // 2D matrix
    int countIslands(vector<vector<int>>);

private:

    int PC(int num);        //periodic conditions on the matrix

    // A function to check if a given
    // cell (row, col) can be included in DFS
    int isSafe(vector<vector<int>> , int row, int col,vector<vector<bool>>);
     
    // A utility function to do DFS for a
    // 2D boolean matrix. It only considers
    // the 8 neighbours as adjacent vertices
    void DFS(vector<vector<int>>, int row, int col, vector<vector<bool>>);

    int ROW,COL;

};

#endif