#include <cmath>
#include <iostream>
#include <fstream>

class YeeGrid {

    public:
        YeeGrid(int Nx, int Ny, double dx, double dy, double dt, double kappa) {
            this->Nx = Nx;
            this->Ny = Ny;
            this->dx = dx;
            this->dy = dy;
            this->dt = dt;
            this->Ex = new double*[Ny];
            this->Ey = new double*[Ny];
            this->Hz = new double*[Ny];
            this->CHx = new double*[Ny];
            this->CHy = new double*[Ny];
            this->CEz = new double*[Ny];
            for (int j = 0; j < Ny; j++) {
                this->Ex[j] = new double[Nx];
                this->Ey[j] = new double[Nx];
                this->Hz[j] = new double[Nx];
                this->CHx[j] = new double[Nx];
                this->CHy[j] = new double[Nx];
                this->CEz[j] = new double[Nx];
            }
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    this->Ex[i][j] = 0;
                    this->Ey[i][j] = 0;
                    this->Hz[i][j] = 0;
                }
            }
            this->kappa = kappa;
            this->tau = dt*10;
            this->t_0 = 6*this->tau;
            this->t = 0;
            this->n = 0;
            this->pulseI = Nx/4;
            this->pulseJ = Ny/4;
        }
        void step();
        void print();
        double getTime();
        int getStep();

    private:
        int Nx;
        int Ny;
        double dx;
        double dy;
        double dt;
        double t;
        int n;
        int pulseI;
        int pulseJ;
        double t_0;
        double tau;
        double kappa;
        double** Ex;
        double** Ey;
        double** Hz;
        double** CEz;
        double** CHx;
        double** CHy;
        void computeCurlEz();
        void computeCurlHx();
        void computeCurlHy();
        void updateEx();
        void updateEy();
        void updateHz();
        void pulse();

};