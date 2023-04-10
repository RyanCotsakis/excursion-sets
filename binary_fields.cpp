#include <Rcpp.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <math.h>       /* sqrt */
using namespace Rcpp;

#define MAP_LENGTH_CHANGE 100


/**
 * Computes map[map[...map[start]...]]
 * 
 * @param map: an array such that map^n[start] is well defined for all n >= 1.
 *  Eventually, map^(n+1)[start] == map^n[start].
 * @param start: the starting index.
 * @return map^(infinity)[start]
 */
int traverse(int map[], int start){
  int next;
  do {
    next = map[start];
    if (next == start){
      return start;
    }
    start = next;
  } while (true);
}


/**
 * Labels connected components (CCs) in a 2D binary field.
 *
 * @param r: The binary field
 * @param connectivity: if 8, CCs connected by corners are considered connected.
 *  If 6, these situations are considered connected with probability 1/2 (stochastic output).
 *  If 4 or else, CCs are not connected by corners.
 * @return a list containing a field like r with the 1s changed to unique labels
 *  that identify the connected component's label, and the same field if 1-r was
 *  passed as an input instead. The pair of rasters respects
 *  the ambiguous connectivity cases. That is, two corners are directly connected
 *  in one raster iff they are not connected in the other raster.
 */
// [[Rcpp::export]]
List cpp_connected_components(NumericMatrix r, int connectivity = 8) {
  int label = 0;
  bool corners_maybe_connected = connectivity == 8 || connectivity == 6;
  bool hole_corners_maybe_connected = connectivity == 4 || connectivity == 6;
  int m = r.ncol(), n = r.nrow();
  int map_length = MAP_LENGTH_CHANGE;
  int *map = new int[map_length];
  map[0] = 0; // empty always maps to empty.
  
  bool components_were_connected;
  bool holes_were_connected;
  
  NumericMatrix ret(n, m);
  NumericMatrix holes(n, m);
  
  for (int row = 0; row < n; row++){
    for (int col = 0; col < m; col++){
      if (r(row,col)){
        
        int neighbours[] = {0, 0, 0, 0};
        if (col > 0){ // We might have a left neighbour
          neighbours[0] = ret(row,col-1);
        }
        if (row > 0){ // We might have a neighbour above
          neighbours[1] = ret(row-1, col);
          if (corners_maybe_connected && col > 0 && (connectivity == 8 || !holes_were_connected)){
            // We might have a neighbour above+left
            neighbours[2] = ret(row-1, col-1);
          }
          if (corners_maybe_connected && col < m-1 && (connectivity == 8 || rand() % 2)){
            // We might have a neighbour above+right
            neighbours[3] = ret(row-1, col+1);
          }
          components_were_connected = bool(neighbours[3]);
        }
        
        int largest_neighbour =  *std::max_element(neighbours, neighbours + 4);
        if (!largest_neighbour){
          // We have no neighbours. Assign a new label
          label++;
          ret(row, col) = label;
          // now store the unique label in the mapping scheme
          if (label == map_length){ // allocate memory for the new label
            map_length += MAP_LENGTH_CHANGE;
            int *temp_map = new int[map_length];
            for (int i = 0; i < map_length - MAP_LENGTH_CHANGE; i++){
              temp_map[i] = map[i];
            }
            delete [] map;
            map = temp_map;
          }
          map[label] = label;
        } else {
          // We have a neighbour.
          ret(row, col) = largest_neighbour;
          int this_neighbour_points_to;
          int max_map = 0;
          for (int i = 0; i < 4; i++){
            this_neighbour_points_to = traverse(map, neighbours[i]);
            if (this_neighbour_points_to > max_map){
              max_map = this_neighbour_points_to;
            }
          }
          for (int i = 0; i < 4; i++){
            if(neighbours[i]){
              map[traverse(map, neighbours[i])] = max_map;
            }
          }
        }
      } else {
        
        int neighbours[] = {0, 0, 0, 0};
        if (col > 0){ // We might have a left neighbour
          neighbours[0] = holes(row,col-1);
        }
        if (row > 0){ // We might have a neighbour above
          neighbours[1] = holes(row-1, col);
          if (hole_corners_maybe_connected && col > 0 && (connectivity == 4 || !components_were_connected)){
            // We might have a neighbour above+left
            neighbours[2] = holes(row-1, col-1);
          }
          if (hole_corners_maybe_connected && col < m-1 && (connectivity == 4 || rand() % 2)){
            // We might have a neighbour above+right
            neighbours[3] = holes(row-1, col+1);
          }
          holes_were_connected = bool(neighbours[3]);
        }
        int largest_neighbour =  *std::max_element(neighbours, neighbours + 4);
        if (!largest_neighbour){
          // We have no neighbours. Assign a new label
          label++;
          holes(row, col) = label;
          // now store the unique label in the mapping scheme
          if (label == map_length){ // allocate memory for the new label
            map_length += MAP_LENGTH_CHANGE;
            int *temp_map = new int[map_length];
            for (int i = 0; i < map_length - MAP_LENGTH_CHANGE; i++){
              temp_map[i] = map[i];
            }
            delete [] map;
            map = temp_map;
          }
          map[label] = label;
        } else {
          // We have a neighbour.
          holes(row, col) = largest_neighbour;
          int this_neighbour_points_to;
          int max_map = 0;
          for (int i = 0; i < 4; i++){
            this_neighbour_points_to = traverse(map, neighbours[i]);
            if (this_neighbour_points_to > max_map){
              max_map = this_neighbour_points_to;
            }
          }
          for (int i = 0; i < 4; i++){
            if(neighbours[i]){
              map[traverse(map, neighbours[i])] = max_map;
            }
          }
        }
      }
    }
  }
  for (int i = 1; i <= label; i++){
    map[i] = traverse(map, i);
  }
  for (int row = 0; row < n; row++){
    for (int col = 0; col < m; col++){
      ret(row, col) = map[int(ret(row, col))];
      holes(row, col) = map[int(holes(row, col))];
    }
  }
  return List::create(Named("cc") = ret,
                      Named("ch") = holes);
}


/**
 * Computes the perimeter of the connected components of a binary field.
 * The boundary of the field is buffered by 0's,
 * and the number of neighbouring pairs who differ (0,1 or 1,0) is returned.
 * 
 * @param r: the binary field
 * @return the number of neighbouring pixel pairs that differ in the buffered field.
 */
// [[Rcpp::export]]
int cpp_perimeter_bierme(NumericMatrix r){
  int count = 0;
  bool previous;
  int m = r.ncol(), n = r.nrow();
  // Scan horizontal changes
  for (int row = 0; row < n; row++){
    previous = r(row, 0);
    for (int col = 1; col < m; col++){
      if (bool(r(row, col)) ^ previous){
        // There has been a change!
        count++;
        previous = !previous;
      }
    }
  }
  // Scan vertical changes
  for (int col = 0; col < m; col++){
    previous = r(0, col);
    for (int row = 1; row < n; row++){
      if (bool(r(row, col)) ^ previous){
        // There has been a change!
        count++;
        previous = !previous;
      }
    }
  }
  return count;
}


/**
 * Computes the perimeter of the connected components of a binary field
 * based on the sum of squared projections.
 * 
 * @param r: the binary field.
 * @param box_size: The pixels are split up into groups of box_size^2 pixels.
 */
// [[Rcpp::export]]
double cpp_perimeter_cotsakis(NumericMatrix r, int box_size){
  int m = r.ncol(), n = r.nrow();
  double total = 0;
  int macro_row = 0;
  int macro_col;
  int previous;
  while (macro_row < n){
    macro_col = 0;
    while (macro_col < m){
      
      // Scan left to right in the box
      int lr_count = 0;
      for (int i = 0; i < box_size && macro_row + i < n; i++){
        if (macro_col == 0){
          previous = r(macro_row+i, 0);
        } else {
          previous = r(macro_row+i, macro_col-1);
        }
        for (int j = 0; j < box_size && macro_col + j < m; j++){
          if (bool(r(macro_row+i, macro_col+j)) ^ previous){
            // There has been a change!
            lr_count++;
            previous = !previous;
          }
        }
      }
      
      // Scan up to down in the box
      int ud_count = 0;
      for (int j = 0; j < box_size && macro_col + j < m; j++){
        if (macro_row == 0){
          previous = r(0, macro_col+j);
        } else {
          previous = r(macro_row-1, macro_col+j);
        }
        for (int i = 0; i < box_size && macro_row + i < n; i++){
          if (bool(r(macro_row+i, macro_col+j)) ^ previous){
            // There has been a change!
            ud_count++;
            previous = !previous;
          }
        }
      }
      total += sqrt(ud_count*ud_count + lr_count*lr_count);
      macro_col += box_size;
    }
    macro_row += box_size;
  }
  return total;
}


/**
 * Computes the "Unbiased" perimeter estimate Given in Appendix B of Bierme et al. 2020.
 * 
 * @param f: the real valued random field.
 * @param u: the threshold at which the excursion set is computed.
 * @return the full perimeter estimate 
 */
// [[Rcpp::export]]
double cpp_perimeter_unbiased(NumericMatrix f, double u){
  int m = f.ncol(), n = f.nrow();
  double total = 0;
  double v[4];
  for (int i = 0; i < n-1; i++){
    for (int j = 0; j < m-1; j++){
      v[0] = f(i, j);
      v[1] = f(i+1, j);
      v[2] = f(i, j+1);
      v[3] = f(i+1, j+1);
      std::sort(v,v+4);
      if (u < v[0] || u > v[3]) continue;
      if (u < v[1]){
        total += (u - v[0])*sqrt(1./pow(v[1]-v[0],2) + 1./pow(v[2]-v[0],2));
      } else if (u < v[2]){
        total += sqrt(1 + pow((u-v[0])/(v[2]-v[0]) - (u-v[1])/(v[3]-v[1]),2));
      } else { //u < v[3]
        total += (v[3] - u)*sqrt(1./pow(v[3]-v[2],2) + 1./pow(v[3]-v[1],2));
      }
    }
  }
  return total;
}


/**
 * Translates an image M over itself and performs a logical operation
 * 
 * @param M: The image
 * @param translations: an (n x 2) matrix of translations: (x,y) or (col,-row)
 * @param operation:
 *    0 -- XOR
 *    1 -- OR
 * @return a new image, same dimensions as the original, but with the operation
 *    applied for each row in translations 
 */
// [[Rcpp::export]]
NumericMatrix cpp_dilation(NumericMatrix M, NumericMatrix translations, int operation){
  int m = M.ncol(), n = M.nrow(), n_trans = translations.nrow();
  NumericMatrix ret(n, m);
  int i1, j1;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      ret(i,j) = 0;
      for (int k = 0; k < n_trans; k++){
        i1 = i - translations(k,1);
        j1 = j + translations(k,0);
        if(0 <= i1 && i1 < n && 0 <= j1 && j1 < m){
          if(operation == 0){ // XOR
            if(M(i,j) != M(i1,j1)){
              ret(i,j) = 1 - ret(i,j);
            }
          } else if(operation == 1){ // OR
            if(M(i,j) || M(i1,j1)){
              ret(i,j) = 1;
              break;
            }
          }
        }
      }
    }
  }
  return ret;
}


/*** R
*/
