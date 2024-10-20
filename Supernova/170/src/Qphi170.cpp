#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib> // 用于随机数
#include <ctime>
#include <fstream>  // 用于文件输出
#include <functional> // 加入这个头文件
#include <random>
#include <chrono> // 包含计时相关的头文件
#include <complex>  // 确保导入复数库

//设置gamma&tau&plambda为常量
const double gamma_phi = 2.0;     // γ，近似为 2
const double tau_phi = 1.0e-2;    // s, 标量粒子的寿命，1 ms
const double lambda_phi = 2.02e6;  // cm, 标量粒子的平均自由程

//单位换算
const double MeV_to_inv_cm = 5.06e10;     // 1 MeV = 5.06e10 cm-1
const double MeV_to_inv_s = 1.52e21;     // 1 MeV = 1.52e21 s-1
const double MeV_to_erg = 1.60e-6;     // 1 MeV = 1.60e-6 erg
const double erg_to_MeV = 6.2415e5; // 1 erg = 6.2415e5 MeV

// 常量定义
const double c = 3.0e10;          // cm/s, 光速
const double hbar = 6.582e-22; // MeV·s, reduced Planck's constant
const double G_F = 1.1663787e-5 * 1e-6; // Fermi coupling constant in MeV^-2
const double PI = 3.141592653589793;    // 圆周率
const double sqrt2 = sqrt(2.0);         // 根号2
const double alpha = 1.0 / 137.0;       // 精细结构常数 α
// 设定初始的 alpha_s 值（在 MZ 质量下）
const double alphas_MZ = 0.118; // Z 玻色子质量下的 αs 值
const double Lambda_QCD = 226.0; // QCD 特征尺度，单位为 MeV

const double m_electron = 0.510998955555;  // 电子 (MeV)
const double m_muon = 105.6583755;         // μ轻子 (MeV)
const double m_tau = 1776.86;              // τ轻子 (MeV)
const double m_up = 2.55;                  // u 夸克 (MeV)
const double m_down = 5.04;                // d 夸克 (MeV)
const double m_strange = 101.0;            // s 夸克 (MeV)
const double m_charm = 1270.0;             // c 夸克 (MeV)
const double m_bottom = 4700.0;            // b 夸克 (MeV)
const double m_top = 173000.0;             // t 夸克 (MeV)

const double MZ = 91187.6; // Z 玻色子质量
const double m_pi = 135.0;        // MeV, pion mass
const double m_N = 938.0;         // MeV, nucleon mass
const double v_Higgs = 246.0 * 1000.0; // 转换为 MeV, Higgs vev

// 输入参数
double f_Tu_p = 0.023;// 胶子项，大约为 0.1
double f_Tu_n = 0.019;// 夸克项，大约为 0.1
double f_Td_p = 0.034;
double f_Td_n = 0.041;
double f_Ts_N = 0.14;

// 超新星相关常量
const double p_F = 200.0;         // MeV, Fermi momentum
const double T_SN = 30.0;         // MeV, supernova temperature
const double R_SN = 10.0 * 1e5;   // 10 km, 转换为 cm

// 函数：计算超新星体积 V_SN
const double V_SN = (4.0 / 3.0) * PI * pow(R_SN, 3);

// 函数：计算 G_phi(u)
double G_phi(double u) {
    double term1 = 1 - (5.0 / 2.0) * pow(u, 2) - (35.0 / 22.0) * pow(u, 4);
    double term2 = (5.0 / 64.0) * (28.0 * pow(u, 3) + 5.0 * pow(u, 5)) * atan(2.0 / u);
    double term3 = (5.0 / 64.0) * sqrt(2.0) * pow(u, 6) / sqrt(2.0 + pow(u, 2)) * atan(2.0 * sqrt(2.0 * (pow(u, 2) + 2)) / pow(u, 2));
    
    return term1 + term2 + term3;
}

// 函数：计算有效耦合到核子 C_phiN
double C_phiN(double sin_theta, double f_Tu_p, double f_Tu_n, double f_Td_p, double f_Td_n, double f_Ts_N) {
    
    // 计算每个夸克的贡献
    double f_Tu_N = f_Tu_p + f_Tu_n; // u 夸克的总贡献
    double f_Td_N = f_Td_p + f_Td_n; // d 夸克的总贡献
    double f_Tq_sum = f_Tu_N + f_Td_N + f_Ts_N; // 夸克的总贡献

    // 计算 f_TG
    double f_TG = 1.0 - f_Tq_sum; // 计算 f_TG

    // 计算有效耦合
    return sqrt(2.0) * sin_theta * m_N / (v_Higgs) * ((6.0 / 27.0) * f_TG + f_Tq_sum); 
}

// 被积函数 - 用于 xi(T, m_phi) 的分子部分
double integrand_numerator(double x, double T, double m_phi) {
    double denominator = exp(x / T) - 1.0;
    if (denominator <= 0) {
        return 1e-10; // 返回一个小正数，防止分母为零
    }

    double sqrt_term = pow(x, 2) - pow(m_phi, 2);
    if (sqrt_term < 0) {
        return 1e-10; // 返回小正数，防止平方根负值
    }

    return x * sqrt(sqrt_term) / denominator;
}

// 被积函数 - 用于 xi(T, m_phi) 的分母部分
double integrand_denominator(double x, double T) {
    double denominator = exp(x / T) - 1.0;
    if (denominator <= 0) {
        return 1e-10; // 返回一个小正数，防止分母为零
    }
    return pow(x, 2) / denominator;
}

// 数值积分 - 使用 Simpson 法来计算积分
double simpson_integral(std::function<double(double, double)> func, double a, double b, int n, double T) {
    double h = (b - a) / n;  // 步长
    double sum = func(a, T) + func(b, T);

    for (int i = 1; i < n; i += 2) {
        sum += 4.0 * func(a + i * h, T);
    }
    for (int i = 2; i < n - 1; i += 2) {
        sum += 2.0 * func(a + i * h, T);
    }

    return (h / 3.0) * sum;
}

// 计算 xi(T, m_phi)
double xi(double T, double m_phi) {
    // 动态调整积分上限，设为 m_phi 的 10 倍
    double upper_limit = 100.0 * m_phi;
    // 积分范围为 [m_phi, upper_limit] 对于分子
    double numerator = simpson_integral(
        [m_phi](double x, double T) { return integrand_numerator(x, T, m_phi); },
        m_phi, upper_limit, 5000, T // 将积分上限设为 upper_limit
    );

    // 积分范围为 [0, upper_limit] 对于分母
    double denominator = simpson_integral(
        integrand_denominator, 
        0.0, upper_limit, 5000, T // 将积分上限设为 upper_limit
    );
    
    // 检查以防止分母为零
    if (denominator == 0) {
        return 0; // 或者根据需要返回一个适当的值
    }

    return numerator / denominator;
}

// 函数：计算逃逸概率 P_esc
double calculate_escape_probability(double R_SN, double gamma_phi, double tau_phi, double lambda_phi) {
    double term1 = exp(-R_SN / (gamma_phi * c * tau_phi)); // exp(-R_SN / γ c τ_φ)
    double term2 = exp(-R_SN / lambda_phi); // exp(-R_SN / λ_φ)
    return term1 * term2;
}


// 对数均匀分布生成随机数的函数
double random_log_uniform(double min, double max, std::default_random_engine& generator) {
    std::uniform_real_distribution<double> distribution(log(min), log(max));
    return exp(distribution(generator));  // 返回对数空间中的随机数并取指数
}


// 主函数
int main() {
    // 全局定义随机数生成器
    std::default_random_engine generator(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
    const int totalIterations = 10000; // 总的迭代次数

    // 开始计时
    auto start = std::chrono::high_resolution_clock::now();

    // 打开 CSV 文件
    std::ofstream outputFile("/home/fyq/DIYbyfyq/Experimentbound/Supernova/170/data/data2.csv");
    if (!outputFile) {
        std::cerr << "无法打开文件!" << std::endl; // 检查文件是否成功打开
        return 1; // 如果失败，返回错误码
    }


    // 设置随机种子
    srand(static_cast<unsigned int>(time(0)));

    // 更新 CSV 文件的表头
    outputFile << "m_phi(GeV),sin_theta^2,Q_phi,C_phiN,G_phi,xi,P_esc,V_SN,lambda_phi,energy loss rate\n";

    std::cout << "开始计算..." << std::endl;

    // 在主函数开头计算 P_esc
    double P_esc = exp(-R_SN / (gamma_phi * c * tau_phi)) * exp(-R_SN / lambda_phi);

    // 计算 G_phi(u)
    double u = m_pi / p_F;
    double G_phi_value = G_phi(u);

    // 生成多个随机点
    for (int i = 0; i < totalIterations; ++i) {
        // 使用对数均匀分布生成 m_phi 和 sin_theta_squared
        double m_phi = random_log_uniform(1.0, 10000.0, generator); // 使用对数均匀分布 m_phi MeV
        double sin_theta_squared = random_log_uniform(1e-13, 1e-3, generator);  // 使用对数均匀分布

        // 从 sin_theta_squared 计算 sin_theta
        double sin_theta = sqrt(sin_theta_squared); // 计算 sin(theta)

        // 计算 C_phiN
        double C_phiN_value = C_phiN(sin_theta, f_Tu_p, f_Tu_n, f_Td_p, f_Td_n, f_Ts_N);

        // 计算 xi(T, m_phi)
        double xi_value = xi(T_SN, m_phi);

        // 计算 Q_phi
        double Q_phi = pow(C_phiN_value, 2) * (11.0 / pow(15.0 * PI, 3)) * pow(T_SN / m_pi, 4) * pow(p_F, 5) * G_phi_value * xi_value;

        // 检查有效性
        if (P_esc < 0 || Q_phi < 0 || lambda_phi < 0 || 
            C_phiN_value < 0 || G_phi_value < 0 || xi_value < 0) {
            continue;  // 跳过此次循环
        }

        // 检查条件 P_esc * Q_phi * V_SN 是否小于等于 10^53
        if (P_esc * Q_phi * V_SN > 0 && P_esc * Q_phi * V_SN * pow(5.06e10, 3) * 1.52e21 * 1.60e-6 > 1e53) {
            double energylossrate = P_esc * Q_phi * V_SN * pow(5.06e10, 3) * 1.52e21 * 1.60e-6; // 计算结果

            // 将满足条件的结果写入 CSV 文件
            outputFile << (m_phi / 1000.0) << ","  // 将 m_phi 从 MeV 转换为 GeV
                       << sin_theta_squared << ","
                       << Q_phi << ","
                       << C_phiN_value << ","
                       << G_phi_value << ","
                       << xi_value << ","
                       << P_esc << ","
                       << V_SN << ","
                       << lambda_phi << ","
                       << energylossrate << "\n";  // 输出 P_esc * Q_phi * V_SN 的乘积
        }

        // 更新进度条
        if (i % (totalIterations / 1000) == 0) { // 每1%更新一次
            std::cout << "\r进度: " << std::fixed << std::setprecision(2)
                      << (static_cast<double>(i) / totalIterations * 100) << "%";
            std::cout.flush();
        }
    }

    // 关闭文件
    outputFile.close();
    std::cout << "\n结果已成功写入到 csv 文件中！" << std::endl;

    // 结束计时
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start; // 计算持续时间
    std::cout << "运行时间: " << duration.count() << " 秒" << std::endl; // 输出运行时间

    return 0; // 确保返回 0
}
