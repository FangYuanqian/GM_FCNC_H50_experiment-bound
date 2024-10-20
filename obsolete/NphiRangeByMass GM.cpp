#include "NphiRangeByMassGM.h" // 确保此文件路径正确
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <sstream>
#include <iomanip>

// 常量定义
const double NPOT = 2e6;       // 总模拟次数，表示模拟中每种类型的介子的总数
const double Nsim = 2e6;       // 总模拟次数，表示所有模拟的总次数
const double d = 128.0;        // 探测器距离目标的距离，单位为米，这里设置为 128 米
const double deltall = 12.0;   // 探测器的长度，单位为米，这里设置为 12 米
const double c = 299792458.0;  // 光速

// 假设P几乎与模型无关，是一个定量（0.5）
const double ξ_K = 1e-2;       // 几何因子
const double P = 0.5;          // 假设的出现概率，表示介子出现的概率，这里设置为 0.5

// 常量定义
const double GF = 1.1663787e-5; 
/// @brief 
/// @brief 
/// @brief 
/// @brief 
const double lambda = 3.1e-7;
/// @brief 
const double A_K_pm = 0.0; 
const double f0_K_pm_pi = 0.96; 
const double m_s = 0.095;  // 单位为 GeV/c^2
const double m_d = 0.005;  // 单位为 GeV/c^2
const double v = 246.0;    // v 的值, 单位：GeV
const double Gamma_K_pm = 0.049; // K^\pm 的衰变宽度，单位为 GeV
const double PI = 3.14159265358979323846;

// 固定值
    double MHH = 10;
    double MH3 = 10;
    double sa = 0.01;
    double M1 = 100;
    double M2 = 0;
    //double MH5 = 10;
    //double sH = 1e-5;
    
//GM
GM mod;
mod.include_partial_two_loop(false);
mod.set_input3(MHH, MH3, MH5, sH, sa, M1, M2);
int K_choice0 = 0;
int K_choice1 = 1;
double xi_H50_W = mod.get_H50_x_W();
std::complex<double> xi_H50_ds = mod.get_H50_xij_down(FIRST, SECOND);
double Br_Kpm= Br_K_PiH50(MH5, xi_H50_W, xi_H50_ds, K_choice0);
double Br_KL= Br_K_PiH50(MH5, xi_H50_W, xi_H50_ds, K_choice1);


// 计算 N_phi 函数
double calculate_N_phi(double M5, double s_H) {
    // 假设的质量和角度值
    double m_K_pm = 0.493;  // K^\pm 质量，单位为 GeV
    double m_pi_pm = 0.140; // π^\pm 质量，单位为 GeV
    double m_KL = 0.497614; 

    // 计算 N_phi
    double N_phi = (NPOT / Nsim) * (N_K_plus * Br_Kpm + N_K_minus * Br_Kpm + N_K_L * Br_KL) * ξ_K * Pphi;

    // 返回 N_phi
    return N_phi;
}

int main() {
    // 示例参数值，实际应根据具体情况设定
    double N_K_plus = NPOT / 3;
    double N_K_minus = NPOT / 3;
    double N_K_L = NPOT / 3;
    double theta_H = 45 * PI / 180;  // 角度转换为弧度

    // 计算 N_phi
    double N_phi = calculate_N_phi();

    // 输出结果
    std::cout << "计算得到的 N_phi: " << N_phi << std::endl;

    return 0;
}
