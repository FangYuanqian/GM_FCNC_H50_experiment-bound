#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <sstream>
#include <iomanip>
#include "Experiment_bound_167.h"
#include <complex>

#include "GM.h"


// 常量定义
const double NPOT = 2e6;       // 总模拟次数，表示模拟中每种类型的介子的总数
const double Nsim = 2e6;       // 总模拟次数，表示所有模拟的总次数
// 假设P几乎与模型无关，是一个定量（0.5）
const double ξ_K = 1e-2;       // 几何因子
const double Pphi = 0.5;          // 假设的出现概率，表示介子出现的概率，这里设置为 0.5
double N_K_plus = 2e6/3;
double N_K_minus = 2e6/3;
double N_K_L = 2e6/3;


// 计算 N_phi 函数
std::complex<double> calculate_N_phi(std::complex<double> Br_Kpm, std::complex<double> Br_KL) { 
    // 计算 N_phi
    std::complex<double> N_phi = (NPOT / Nsim) * (N_K_plus * Br_Kpm + N_K_minus * Br_Kpm + N_K_L * Br_KL) * ξ_K * Pphi;
    return N_phi;
}
