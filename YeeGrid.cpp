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

void YeeGrid::step() {
    computeCurlEz();
    updateHz();
    pulse();
    computeCurlHx();
    computeCurlHy();
    updateEx();
    updateEy();
    n++;
    t += dt;
}

double YeeGrid::getTime() {
    return t;
}

int YeeGrid::getStep() {
    return n;
}

void YeeGrid::pulse() {
    Hz[pulseI][pulseJ] = Hz[pulseI][pulseJ] + std::exp(-pow((t-t_0)/tau,2));
}

void YeeGrid::print() {
    std::ofstream exFile, eyFile, hzFile;
    std::string nxn = std::to_string(Nx) + "x" + std::to_string(Ny);
    std::string path = "results/" + nxn + "/";
    // std::cout << path << "\n";
    exFile.open(path + "Ex_" + std::to_string(n) + ".csv");
    eyFile.open(path + "Ey_" + std::to_string(n) + ".csv");
    hzFile.open(path + "Hz_" + std::to_string(n) + ".csv");
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            exFile << std::to_string(Ex[i][j]) + ",";
            eyFile << std::to_string(Ey[i][j]) + ",";
            hzFile << std::to_string(Hz[i][j]) + ",";
        }
        exFile << "\n";
        eyFile << "\n";
        hzFile << "\n";
    }
    exFile.close();
    eyFile.close();
    hzFile.close();
}

void YeeGrid::updateEx() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            Ex[i][j] = Ex[i][j] + dt*pow(kappa,2)*(CHx[i][j]);
        }
    }
}

void YeeGrid::updateEy() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            Ey[i][j] = Ey[i][j] + dt*pow(kappa,2)*(CHy[i][j]);
        }
    }
}

void YeeGrid::updateHz() {
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            Hz[i][j] = Hz[i][j] + -dt*(CEz[i][j]);
        }
    }
}

void YeeGrid::computeCurlEz() {
    for (int i = 0; i < Nx-1; i++) {
        for (int j = 0; j < Ny-1; j++) {
            CEz[i][j] = (Ey[i+1][j] - Ey[i][j])/dx - (Ex[i][j+1] - Ex[i][j])/dy;
        }
        CEz[i][Ny-1] = (Ey[i+1][Ny-1] - Ey[i][Ny-1])/dx - (Ex[i][0] - Ex[i][Ny-1])/dy;
    }
    for (int j = 0; j < Ny-1; j++) {
        CEz[Nx-1][j] = (Ey[0][j] - Ey[Nx-1][j])/dx - (Ex[Nx-1][j+1] - Ex[Nx-1][j])/dy;
    }
    CEz[Nx-1][Ny-1] = (Ey[0][Ny-1] - Ey[Nx-1][Ny-1])/dx - (Ex[Nx-1][0] - Ex[Nx-1][Ny-1])/dy;
}

void YeeGrid::computeCurlHx() {
    for (int i = 0; i < Nx; i++) {
        CHx[i][0] = (Hz[i][0] - Hz[i][Ny-1])/dy;
        for (int j = 1; j < Ny; j++) {
            CHx[i][j] = (Hz[i][j] - Hz[i][j-1])/dy;
        }
    }
}

void YeeGrid::computeCurlHy() {
    for (int j = 0; j < Ny; j++) {
        CHy[0][j] = (Hz[0][j] - Hz[Nx-1][j])/dy;
        for (int i = 1; i < Nx; i++) {
            CHy[i][j] = -(Hz[i][j] - Hz[i-1][j])/dy;
        }
    }
}