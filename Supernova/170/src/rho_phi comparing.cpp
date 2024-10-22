#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib> 
#include <ctime>   
#include <fstream>  
#include <functional> 
#include <random>
#include <chrono> 
#include <complex>  
#include <vector>
#include <string>
#include <limits>

// 单位换算
const double MeV_to_inv_cm = 5.06e10;     // 1 MeV = 5.06e10 cm-1
const double MeV_to_inv_s = 1.52e21;     // 1 MeV = 1.52e21 s-1
const double MeV_to_erg = 1.60e-6;     // 1 MeV = 1.60e-6 erg
const double erg_to_MeV = 6.2415e5;    // 1 erg = 6.2415e5 MeV
const double PI = 3.141592653589793;   // 圆周率

// 超新星相关常量
const double T_SN = 30.0;   // MeV, 超新星温度
const double R_SN = 10.0 * 1e5;   // 10 km, 转换为 cm

// 被积函数
double integrand(double E, double m_phi) {
    double result = (sqrt(E * E - m_phi * m_phi) / (exp(E / T_SN) + 1)) * E * E;
    return result;
}

// 能量密度计算
double calculate_energy_density(double m_phi) {
    double integral = 0.0;
    const int n = 1000; // 积分分区数量
    const double E_max = 10 * T_SN; // 设定合适的上限
    double h = (E_max - m_phi) / n; // 使用 m_phi 作为下限

    for (int i = 0; i < n; i++) {
        double E = m_phi + i * h; // 使用 m_phi 作为起始点
        double integrand_value = integrand(E, m_phi);
        integral += integrand_value;
    }

    integral *= h; // 积分结果乘以步长

    // 乘以前面的常数因子
    return (1 / (2 * PI * PI)) * integral; // 使用 g 计算最终结果
}

// 使用公式计算 rho_phi
double calculate_rho_phi(double T_SN) {
    return (pow(PI, 2) * pow(T_SN, 4)) / 30.0;
}

int main() {
    // 直接输入数据
    std::vector<double> m_phi_values = {0.001, 0.01, 0.1}; 
    std::vector<double> sin_theta_squared_values = {1e-12, 1e-11, 1e-10, 1e-9};

    for (size_t i = 0; i < m_phi_values.size(); ++i) {
        double m_phi = m_phi_values[i] * 1e3; // 将 GeV 转换为 MeV
        double sin_theta_squared = sin_theta_squared_values[i];

        std::cout << "m_phi: " << m_phi_values[i] << " GeV, sin_theta_squared: " << sin_theta_squared 
                  << std::endl;

        // 计算积分方法的能量密度
        double rho_phi_integral = calculate_energy_density(m_phi);

        // 使用公式直接计算
        double rho_phi_formula = calculate_rho_phi(T_SN);

        // 输出对比
        std::cout << "积分法 rho_phi: " << std::setprecision(6) << std::fixed << rho_phi_integral 
                  << " MeV^4, 公式法 rho_phi: " << std::setprecision(6) << std::fixed 
                  << rho_phi_formula << " MeV^4" << std::endl;
    }

    return 0;
}
