#pragma once
#include<iostream>
#include"LinearEquation.h"
#include <Eigen/Dense>
#include <chrono>  // 用于测量时间
#include <random>  // 用于生成随机数

using namespace std;
using namespace Eigen;

// 生成随机矩阵和向量
void generate_random_data(MatrixXd &A, VectorXd &b) {
    random_device rd;  // 获取随机数种子
    mt19937 gen(rd()); // 初始化随机数生成器
    uniform_real_distribution<> dis(-10.0, 10.0);  // 生成范围为 [-10, 10] 的随机数

    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            A(i, j) = dis(gen);  // 生成随机系数
        }
        b(i) = dis(gen);  // 生成随机右边的向量
    }
}

int main() {
    // 随机生成矩阵和向量的规模
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> size_dis(2, 10);  // 矩阵和向量的规模范围为 [2, 10]

    int n = size_dis(gen);  // 随机规模
    MatrixXd A(n, n);
    VectorXd b(n);

    // 生成随机的矩阵和向量
    generate_random_data(A, b);

    // 测试自定义的 LinearEquation 类
    LinearEquation<double> lq(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            lq.setMatrix(A(i, j), i, j);
        }
    }
    for (int i = 0; i < n; i++) {
        lq.setResult(b(i), i);
    }

    cout << "Matrix A:\n" << A << endl;
    cout << "Vector b:\n" << b << endl;

    cout << "Before elimination: " << endl;
    lq.show();

    auto start_lq = chrono::high_resolution_clock::now();
    lq.gaussian_elimination();


    cout << "After elimination: " << endl;
    vector<double> result_lq = lq.solution();
    auto end_lq = chrono::high_resolution_clock::now();
    auto duration_lq = chrono::duration_cast<chrono::microseconds>(end_lq - start_lq).count();
    for (int i = 0; i < result_lq.size(); i++) {
        cout << result_lq[i] << " ";
    }
    cout << endl;
    cout << "Custom LinearEquation class took " << duration_lq << " microseconds." << endl;

    // 使用 Eigen 库求解（使用矩阵逆）
    auto start_eigen = chrono::high_resolution_clock::now();
    VectorXd x = A.inverse() * b;  // 使用矩阵逆求解
    auto end_eigen = chrono::high_resolution_clock::now();
    auto duration_eigen = chrono::duration_cast<chrono::microseconds>(end_eigen - start_eigen).count();

    cout << "Eigen solution: " << x.transpose() << endl;
    cout << "Eigen took " << duration_eigen << " microseconds." << endl;

    return 0;
}
