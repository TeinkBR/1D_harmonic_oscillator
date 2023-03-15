#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

class Emulation {
public:
    Emulation(const std::vector<Eigen::MatrixXd> &trainingList) : trainingList(trainingList) {
        // Initialize hbar, mass, e_max, meshpoint, weight, basis_wave_function_ilevel_ix, and energy_level
        // based on the values from the first item in the training list.
        // Please replace the placeholders with the actual values.
        
        hbar = 1.0; // Replace with the actual value
        mass = 1.0; // Replace with the actual value
        e_max = 1.0; // Replace with the actual value
        meshpoint = Eigen::VectorXd::LinSpaced(100, -10, 10); // Replace with the actual value
        weight = 1.0; // Replace with the actual value
        basis_wave_function_ilevel_ix = Eigen::MatrixXd::Random(10, 10); // Replace with the actual value
        energy_level = Eigen::VectorXd::LinSpaced(10, 0, 9); // Replace with the actual value

        computeNMatrix();
        computeHMatrix();
    }

    void startEmu(double emulationOmega, const std::function<double(double)> &emulationVsFuncOfXUnscaled) {
        // Compute the scaling factor and V_s
        double scalingFactorB = hbar / sqrt(mass * emulationOmega);
        Eigen::VectorXd V_s_ix = meshpoint.unaryExpr([=](double x) { return emulationVsFuncOfXUnscaled(x * scalingFactorB) / emulationOmega; });

        // Compute V_s_ilevel_ilevel
        Eigen::MatrixXd V_s_ilevel_ilevel = basis_wave_function_ilevel_ix * V_s_ix.asDiagonal() * basis_wave_function_ilevel_ix.transpose();

        // Compute emulation_H
        Eigen::MatrixXd emulation_H = H_p_sqr_r_sqr_ilevel_ilevel + V_s_ilevel_ilevel;

        // Compute H_matrix
        H_matrix.resize(trainingList.size(), std::vector<Eigen::MatrixXd>(trainingList.size()));
        for (size_t i = 0; i < trainingList.size(); ++i) {
            for (size_t j = 0; j < trainingList.size(); ++j) {
                H_matrix[i][j] = trainingList[i].col(0).conjugate() * emulation_H * trainingList[j].col(0);
            }
        }

        // Solve the generalized eigenvalue problem
        Eigen::MatrixXd H_mat = toEigenMatrix(H_matrix);
        Eigen::MatrixXd N_mat = toEigenMatrix(N_matrix);
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_mat, N_mat);

        // Get the eigenvalues and eigenvectors
        eigvalEmuH = eigensolver.eigenvalues();
        eigvecEmuH = eigensolver.eigenvectors();
    }

    Eigen::VectorXd getEigenvalues() const { return eigvalEmuH; }
    Eigen::MatrixXd getEigenvectors() const { return eigvecEmuH; }

private:
    void computeNMatrix() {
        // Compute N_matrix based on the training list
        // Please replace the placeholders with the actual values.
    }

    void computeHMatrix() {
       
