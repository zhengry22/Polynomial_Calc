#pragma once
#include<iostream>
#include <cassert>
using namespace std;

template<typename T>
class LinearEquation {
    /*
        This class is used to solve the linear equation Ax = b.
    */
private:
    size_t rows, cols;
    T** matrix; // A
    T* result; // b
    inline void swap(size_t a, size_t b) {
        // swap row
        T* tmp = this->matrix[a];
        this->matrix[a] = this->matrix[b];
        this->matrix[b] = tmp;
        // swap vec
        T tmp2 = this->result[a];
        this->result[a] = this->result[b];
        this->result[b] = tmp2;
    }
    inline void to_1(size_t x, size_t y) {
        /*
            reduct the x-th row, so that the y-th col
            of the x-th col could be 1.
        */
        assert(x < this->rows && "Error: x to large in to_1!");
        assert(this->matrix[x][y] != 0 && "Error: pivot is 0!");
        for (int i = 0; i < this->cols; i++) {
            this->matrix[x][i] /= this->matrix[x][y];
        }
        this->result[x] /= this->matrix[x][y];
    }
    inline void row_reduct(size_t x, size_t y, T factor) {
        /*
            the whole x-th row needs to minus difference.
        */
        assert(x < this->rows && "Error: x to large in to_1!");
        assert(y < this->rows && "Error: y to large in to_1!");
        for (int i = 0; i < this->cols; i++) {
            this->matrix[x][i] -= this->matrix[y][i] * factor;
        }
        this->result[x] -= this->result[y] * factor;
    }
public:
    LinearEquation(size_t x, size_t y):rows(x),cols(y){
        this->matrix = new T*[x];
        for (int i = 0; i < x; i++) {
            this->matrix[i] = new T[y];
        }
        this->result = new T[x];
    }
    ~LinearEquation() {
        for (int i = 0; i < this->rows; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;
        delete[] result;
    }
    inline bool isSquare() {return this->rows == this->cols;}
    inline void setMatrix(T coeff, size_t x, size_t y) {
        assert(x < this->rows && "Error: exceed maximum row");
        assert(y < this->cols && "Error: exceed maximum column");
        this->matrix[x][y] = coeff;
    }
    inline void setResult(T coeff, size_t x) {
        assert(x < this->rows && "Error: exceed maximum row");
        this->result[x] = coeff;
    }
    inline void show() {
        cout << "A: " << endl;
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->cols; j++) {
                cout << this->matrix[i][j] << " ";
            }
            cout << endl;
        }
        cout << "B: " << endl;
        for (int i = 0; i < this->rows; i++) {
            cout << this->result[i] << " ";
        } 
        cout << endl;
    }
    void gaussian_elimination() {
        int row_id = 0;
        for (int i = 0; i < this->rows; i++) {
            /*
                Note that the i-th col of the i-th row maybe 0 after elimination.
            */
            if (matrix[row_id][i] == 0) {
                // swap
                bool swap_flag = false;
                for (int j = row_id; j < this->rows; j++) {
                    if (matrix[j][i] != 0) {
                        // swap row i and j for the matrix and vector
                        swap(i, j);
                        swap_flag = true;
                        break;
                    }
                }
                if (!swap_flag) {
                    /*
                        This indicates that the i-th col is all zeros,
                        and therefore can not be swapped. So we move on
                        to the next row and col.
                    */
                    continue;
                }
            } 
            /*
                Now the first row should be non zero. First we eliminate the 
                pivot to 1 and then eliminate.
            */
            to_1(row_id, i);
            for (int j = 0; j < this->rows; j++) {
                if (j == row_id) continue;
                T factor = this->matrix[j][i] / this->matrix[row_id][i];
                row_reduct(j, row_id, factor);
            }
            row_id++;
        }
    }
};