#include "/home/fyq/H5/include/NphiRangeByMassGM.h" // 确保此文件路径正确
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

// 计算 \xi_{\phi}^{ds}
double calculate_xi_phi_ds(double m_H5, double theta_H) {
    // 假设值
    return 1.0;  
}

// 计算 \xi_W^5
double calculate_xi_W(double m_H5, double theta_H) {
    // 假设值
    return 1.0; 
}

// 计算分支比 Br(K^\pm \rightarrow \pi^\pm H_5^0)
double calculate_Br_K_pm(double m_K_pm, double m_pi_pm, double m_H5, double theta_H) {
    // 假设值
    return 1.0;  // 所有的假设值都设为 1
}

// 计算分支比 Br(K_L \to \pi^0 H_5^0)
double calculate_Br_K_L(double m_K_L, double m_pi_0, double m_H5, double theta_H) {
    // 假设值
    return 1.0;  // 所有的假设值都设为 1
}

// 计算 N_phi 函数
double calculate_N_phi(double N_K_plus, double N_K_minus, double N_K_L, double theta_H) {
    // 假设的质量和角度值
    double m_K_pm = 0.493;  // K^\pm 质量，单位为 GeV
    double m_pi_pm = 0.140; // π^\pm 质量，单位为 GeV
    double m_H5 = -0.4562248628884826; // 选择的 m_H5 值,范围：(-0.3675694698354661, -0.5448802559414991)

    // 计算每种介子的衰变分支比
    double Br_K_plus = calculate_Br_K_pm(m_K_pm, m_pi_pm, m_H5, theta_H);
    double Br_K_minus = calculate_Br_K_pm(m_K_pm, m_pi_pm, m_H5, theta_H);
    double Br_K_L = calculate_Br_K_L(m_K_pm, m_pi_pm, m_H5, theta_H);

    // 计算 N_phi
    double N_phi = (NPOT / Nsim) * (N_K_plus * Br_K_plus + N_K_minus * Br_K_minus + N_K_L * Br_K_L) * ξ_K * P;

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
    double N_phi = calculate_N_phi(N_K_plus, N_K_minus, N_K_L, theta_H);

    // 输出结果
    std::cout << "计算得到的 N_phi: " << N_phi << std::endl;

    return 0;
}
