//CONTEGGIO ISOLE IN UNA MATRICE
//23MAGGIO2021

//g++ -std=c++11 Island.cpp -o Island.x


#include <iostream>
#include <fstream>

#define ROW 64
#define COL 64

using namespace std;

int PC(int num){           
  if (num == ROW){num = 0;}
  else if (num == -1){num = ROW-1;}
  return num;
}

// A function to check if a given
// cell (row, col) can be included in DFS
int isSafe(int M[][COL], int row, int col, bool visited[][COL])
{
    // row number is in range, column
    // number is in range and value is 1
    // and not yet visited
    return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL) && (M[row][col] && !visited[row][col]);

}
 
// A utility function to do DFS for a
// 2D boolean matrix. It only considers
// the 8 neighbours as adjacent vertices
void DFS(int M[][COL], int row, int col, bool visited[][COL])
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
        if (isSafe(M, PC(row + rowNbr[k]), PC(col + colNbr[k]), visited)){
            DFS(M, PC(row + rowNbr[k]), PC(col + colNbr[k]), visited);
        }
}
 
// The main function that returns
// count of islands in a given boolean
// 2D matrix
int countIslands(int M[][COL]){

    // Make a bool array to mark visited cells.
    // Initially all cells are unvisited
    bool visited[ROW][COL];
    memset(visited, 0, sizeof(visited));
 
    // Initialize count as 0 and
    // travese through the all cells of
    // given matrix
    int count = 0;

    for (int i = 0; i < ROW; ++i){
        for (int j = 0; j < COL; ++j){
 
            // If a cell with value 1 is not  
            // visited yet, then new island found
            // Visit all cells in this island.
            if (M[i][j] && !visited[i][j]) {
               
                DFS(M, i, j, visited);
 
                // and increment island count
                ++count;
            }
        }
    }    
    return count;
}


int main(){

    ifstream filein("matriceFin.dat");
    int M[ROW][COL];
    int appo;
    for(int i=0;i<ROW;++i){
        for(int j=0;j<COL;j++){
            filein>>M[i][j];
        }
    }


    cout << "nIsole=" << countIslands(M)<<endl;



    return 0;
}