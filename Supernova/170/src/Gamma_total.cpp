#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

const double G_F = 1.1663787e-5 * 1e-6; // Fermi coupling constant in MeV^-2
const double PI = 3.141592653589793;    // 圆周率
const double sqrt2 = sqrt(2.0);         // 根号2
const double alpha = 1.0 / 137.0;       // 精细结构常数 α
const double m_pi = 0.13957 * 1e3;      // π 粒子的质量，单位是 MeV
const double m_e = 0.000511 * 1e3;      // 电子的质量，单位为 MeV
const double m_mu = 0.10566 * 1e3;      // μ子质量，单位为 MeV

// 定义 f(tau) 函数
std::complex<double> f(double tau) {
    if (tau <= 1) {
        return std::pow(asin(sqrt(tau)), 2); // tau ≤ 1 时，计算反正弦的平方
    } else {
        double term1 = log((1 + sqrt(1 - 1/tau)) / (1 - sqrt(1 - 1/tau))); // 计算 log
        return -0.25 * std::pow(term1, 2) - std::complex<double>(0, 0.25 * PI * PI); // 返回复数结果
    }
}

// 定义 A_1/2(tau) 函数
double A_1_2(double tau) {
    return 2.0 * (tau + (tau - 1) * f(tau).real()) / (tau * tau); // 计算 A_1/2
}

// 定义 A_1(tau) 函数
double A_1(double tau) {
    return - (2 * tau * tau + 3 * tau + 3 * (2 * tau - 1) * f(tau).real()) / (tau * tau); // 计算 A_1
}

// 计算衰变宽度 Γ(φ → f \bar{f})
double calculate_Gamma_phi_to_f(double m_phi, double m_f, double s_theta) {
    double term1 = sqrt(1 - (4 * pow(m_f, 2) / pow(m_phi, 2))); // 计算动量因子
    return s_theta * s_theta * (G_F * pow(m_f, 2) * m_phi) / (4 * sqrt2) * pow(term1, 1.5); // 计算衰变宽度
}

// 计算衰变宽度 Γ(φ → γγ)
double calculate_Gamma_phi_to_gamma(double m_phi, double s_theta) {
    // 计算常数项
    double term1 = s_theta * s_theta * G_F * pow(alpha, 2) * pow(m_phi, 3) / (128 * sqrt2 * pow(PI, 3));
    double sum = 0.0; // 初始化总和

    // 定义费米子数量、颜色数和电荷等
    std::vector<double> N_c = {1, 3, 3};                     // 电子、u 夸克、d 夸克的颜色数
    std::vector<double> Q_f = {1.0, 2.0 / 3.0, -1.0 / 3.0};  // 电子、u 夸克、d 夸克的电荷
    std::vector<double> m_f = {0.000511 * 1e3, 0.002 * 1e3, 0.004 * 1e3}; // 电子、u 夸克、d 夸克的质量 (GeV)

    // 遍历所有费米子
    for (size_t i = 0; i < N_c.size(); ++i) {
        // 计算 tau
        double tau_f = m_phi * m_phi / (4 * pow(m_f[i], 2)); // 用质量而不是电荷计算 tau
        sum += N_c[i] * pow(Q_f[i], 2) * A_1_2(tau_f);       // 累加 A_1/2 的贡献
    }

    // 加入 A_1(tau)
    // 避免除以零错误
    double tau = (m_phi != 0) ? m_phi * m_phi / (4 * 1e-12) : 0.0; // 假设光子质量接近零但不为零
    sum += A_1(tau); // 加入 A_1 的贡献

    // 计算并返回总衰变宽度
    return term1 * pow(sum, 2);
}

// 强耦合常数 alpha_s 的定义（这里你需要根据 m_phi 实现这个函数）
double alpha_s(double m_phi) {
    // 在某个 m_phi 下的 alpha_s 值。此处你可能需要引入 QCD 运行公式进行更精确的计算。
    return 0.3; // 暂时设定为常数
}

// 计算衰变宽度 Γ(φ → gg)
double calculate_Gamma_phi_to_glue(double m_phi, double s_theta) {
    double alpha_s_value = alpha_s(m_phi); // 获取强耦合常数
    double sum = 0.0;                      // 初始化总和

    // 假设考虑的夸克质量（单位：MeV）
    std::vector<double> quark_masses = {0.002 * 1e3, 0.004 * 1e3, 0.095 * 1e3}; // u, d, s 夸克的质量
    std::vector<double> N_c = {3.0, 3.0, 3.0};                                  // 每种夸克的颜色数

    // 累加每种夸克的贡献
    for (size_t i = 0; i < quark_masses.size(); ++i) {
        double m_q = quark_masses[i];
        double tau_q = m_phi * m_phi / (4 * m_q * m_q); // 计算 tau
        sum += N_c[i] * A_1_2(tau_q);                  // 累加 A_1/2 的贡献
    }

    // 计算并返回总衰变宽度
    double width = s_theta * s_theta * G_F * pow(alpha_s_value, 2) * pow(m_phi, 3) / (36 * sqrt2 * pow(PI, 3));
    width *= pow((3.0 / 4.0 * sum), 2); // 计算最终的宽度

    return width;
}

// 计算总衰变宽度
double calculate_Gamma_total(double m_phi, double s_theta) {
    double Gamma_total = 0.0;

    // 计算光子衰变
    Gamma_total += calculate_Gamma_phi_to_gamma(m_phi, s_theta);

    // 如果 m_phi > 2m_π，计算所有衰变
    if (m_phi > 2 * m_pi) {
        Gamma_total += calculate_Gamma_phi_to_f(m_phi, m_e, s_theta); // e+e- 渠道
        Gamma_total += calculate_Gamma_phi_to_f(m_phi, m_mu, s_theta); // μ+μ- 渠道
        Gamma_total += calculate_Gamma_phi_to_glue(m_phi, s_theta);   // 胶子衰变
    } else {
        // 否则，只计算 e+e- 和 μ+μ- 渠道
        if (m_phi > 2 * m_e) {
            Gamma_total += calculate_Gamma_phi_to_f(m_phi, m_e, s_theta); // e+e- 渠道
        }
        if (m_phi > 2 * m_mu) {
            Gamma_total += calculate_Gamma_phi_to_f(m_phi, m_mu, s_theta); // μ+μ- 渠道
        }
    }

    return Gamma_total;
}

int main() {
    double m_phi = 1.0;   // φ 粒子的质量 (单位：GeV)
    double s_theta = 0.1; // 假设的混合角

    // 计算并输出总衰变宽度
    double Gamma_total = calculate_Gamma_total(m_phi, s_theta);
    std::cout << "φ 粒子的总衰变宽度为：" << Gamma_total << " MeV" << std::endl;

    return 0;
}
