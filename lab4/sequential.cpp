#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <algorithm>

const double a = 1e5; // параметр уравнения
const double epsilon = 1e-8; // порог сходимости
const double D_x = 2.0, D_y = 2.0, D_z = 2.0; // область моделирования: [-1, 1] x [-1, 1] x [-1, 1]
const double x_0 = -1.0, y_0 = -1.0, z_0 = -1.0;
const int Nx = 100, Ny = 100, Nz = 100; // параметры сетки
const double hx = D_x / (Nx - 1); // шаги сетки по различным осям
const double hy = D_y / (Ny - 1);
const double hz = D_z / (Nz - 1);
const int max_iterations = 10000;

inline double phi(double x, double y, double z) {
    return x*x + y*y + z*z;
}

inline double rho(double x, double y, double z) {
    return 6.0 - a * phi(x, y, z);
}

int idx(int i, int j, int k) {
    return i * Ny * Nz + j * Nz + k;
}

int main() {
    std::vector<double> phi_old(Nx * Ny * Nz, 0.0);
    std::vector<double> phi_new(Nx * Ny * Nz, 0.0);
    std::vector<double> rho_vals(Nx * Ny * Nz, 0.0);

    // граничные условия и начальное приближение
    for (int i = 0; i < Nx; ++i) {
        double x = x_0 + i * hx;
        for (int j = 0; j < Ny; ++j) {
            double y = y_0 + j * hy;
            for (int k = 0; k < Nz; ++k) {
                double z = z_0 + k * hz;
                rho_vals[idx(i,j,k)] = rho(x, y, z);
                // граница: i=0, i=Nx-1, j=0, j=Ny-1, k=0, k=Nz-1
                if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1 || k == 0 || k == Nz-1) {
                    phi_old[idx(i,j,k)] = phi(x,y,z);
                }
            }
        }
    }

    int iterations = 0;
    double max_diff = 0.0;

    auto start = std::chrono::high_resolution_clock::now();

    // итерационный процесс якоби по внутренним узлам
    do {
        max_diff = 0.0;
        for (int i = 1; i < Nx-1; ++i) {
            for (int j = 1; j < Ny-1; ++j) {
                for (int k = 1; k < Nz-1; ++k) {
                    // формула для итерационного процесса метода Якоби
                    double term_x = (phi_old[idx(i+1,j,k)] + phi_old[idx(i-1,j,k)]) / (hx*hx);
                    double term_y = (phi_old[idx(i,j+1,k)] + phi_old[idx(i,j-1,k)]) / (hy*hy);
                    double term_z = (phi_old[idx(i,j,k+1)] + phi_old[idx(i,j,k-1)]) / (hz*hz);

                    // координаты точки в физическом пространстве - для вычисления фи
                    double x = x_0 + i * hx;
                    double y = y_0 + j * hy;
                    double z = z_0 + k * hz;

                    // формула якоби
                    double numerator = term_x + term_y + term_z - rho_vals[idx(i,j,k)];
                    double denom = 2.0/(hx*hx) + 2.0/(hy*hy) + 2.0/(hz*hz) + a;
                    phi_new[idx(i,j,k)] = numerator / denom;

                    double diff = std::fabs(phi_new[idx(i,j,k)] - phi_old[idx(i,j,k)]);
                    max_diff = std::max(max_diff, diff);
                }
            }
        }

        // обновляем старое приближение
        for (int i = 1; i < Nx-1; ++i) {
            for (int j = 1; j < Ny-1; ++j) {
                for (int k = 1; k < Nz-1; ++k) {
                    phi_old[idx(i,j,k)] = phi_new[idx(i,j,k)];
                }
            }
        }

        ++iterations;
    } while (max_diff > epsilon && iterations < max_iterations);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "итерации: " << iterations << std::endl;
    std::cout << "max diff: " << max_diff << "; epsilon: " << epsilon << std::endl;
    std::cout << "time: " << elapsed.count() << std::endl;

    return 0;
}

/*
компилируй так:
    g++ -O3 -o sequential sequential.cpp
    ./sequential
*/
