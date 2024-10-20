#include "GM.h"
#include "Experiment_bound_167.h"
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "Constants.h"
#include "FCNC.h"
#include "LTLauncher.h"
#include "SM.h"
#include "forward_scalar_rate.h"
#include "lorentz_vector.h"
#include "meson_decay.h"
#include "meson_rate.h"
#include "spdlog_wrapper.h"
#include "utilities.h"
#include <cmath>
#include <complex>


int main() {
    // 常量定义
    double MHH = 10;
    double MH3 = 10;
    double sa = 0.01;
    double M1 = 100;
    double M2 = 0;
    double MH5 = 10;
    double sH = 1e-5;

    // 打开文件
    std::ofstream outFile("/home/wang/H5_Gamma/test_167.csv", std::ios::app);
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // 如果文件是空的，写入表头
    if (outFile.tellp() == 0) {
        outFile << "MH5, sH, N_phi_GM" << std::endl;
    }

    // 假设 GM 是一个已定义的类
    GM mod;
    mod.include_partial_two_loop(false);
    mod.set_input3(MHH, MH3, MH5, sH, sa, M1, M2);
    int K_choice0 = 0;
    int K_choice1 = 1;
    double xi_H50_W = mod.get_H50_x_W();
    std::complex<double> xi_H50_ds = mod.get_H50_xij_down(FIRST, SECOND);
    std::complex<double> Br_Kpm = Br_K_PiH50(MH5, xi_H50_W, xi_H50_ds, K_choice0);
    std::complex<double> Br_KL = Br_K_PiH50(MH5, xi_H50_W, xi_H50_ds, K_choice1);
    std::complex<double> N_phi_GM = calculate_N_phi(Br_Kpm, Br_KL);

    outFile << MH5 << ","
            << sH << ","
            << N_phi_GM << std::endl;

    std::cout << "Data written to /home/wang/H5_Gamma/test_167.csv" << std::endl;

    outFile.close();
    return 0;
}