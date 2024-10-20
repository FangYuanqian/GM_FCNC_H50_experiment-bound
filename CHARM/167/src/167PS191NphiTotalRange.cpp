#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath> // 引入数学库

// 常量定义
const double NPOT = 2e6;       // 总模拟次数，表示模拟中每种类型的介子的总数
const double Nsim = 2e6;       // 总模拟次数，表示所有模拟的总次数
const double d = 128.0;        // 探测器距离目标的距离，单位为米，这里设置为 128 米
const double deltall = 12.0;   // 探测器的长度，单位为米，这里设置为 12 米
const double c = 299792458.0;  // 光速，单位为米每秒

// 假设Nk&ξk比例系数（与实验无关），取为定值
const double ξ_K = 1e-2;       // 几何因子，用于计算某些物理过程中的几何效应
const double Pphi = 0.5;       // 假设的出现概率，表示介子出现的概率

// 计算标量粒子在K介子静止系中的动量
double p_phi_CM(double m_phi, double m_pi, double m_K) {
    double term1 = (m_phi + m_pi) / m_K;
    double term2 = (m_phi - m_pi) / m_K;
    
    double arg1 = 1 - pow(term1, 2);
    double arg2 = 1 - pow(term2, 2);

    // 如果arg1或arg2小于0，返回0
    if (arg1 < 0 || arg2 < 0) {
        return 0;
    }

    return (m_K / 2) * sqrt(arg1 * arg2);
}

// 计算每种介子的衰变分支比 Br(K → πφ)
double Br_K_to_πφ(const std::string& K_type, double θ, double p_phi_CM_value) {
    // 定义不同类型 K 介子的质量
    double m_K;
    if (K_type == "K+") {
        m_K = 0.493677;  // K+ 的质量，单位：GeV/c^2
    } else if (K_type == "K-") {
        m_K = 0.493677;  // K- 的质量，单位：GeV/c^2
    } else if (K_type == "K_L") {
        m_K = 0.497611;  // K_L 的质量，单位：GeV/c^2
    } else {
        return 0; // 如果介子类型不匹配，返回 0
    }

    // 根据介子类型的质量计算分支比
    if (K_type == "K_L") {
        return 5.7e-3 * (2 * p_phi_CM_value / m_K) * pow(θ, 2);
    } else {
        return 1.6e-3 * (2 * p_phi_CM_value / m_K) * pow(θ, 2);
    }
}

// 计算 N_phi
double calculate_N_phi(double m_phi, double θ_squared) {
    double θ = std::sqrt(θ_squared);  // 计算 θ
    double m_pi = 0.13957018;          // π 的质量，单位：GeV/c^2
    double p_phi_CM_value = p_phi_CM(m_phi, m_pi, 0.493677); // 使用 K+ 的质量作为示例
    double N_K_plus = NPOT / 3;       // K+ 的数量
    double N_K_minus = NPOT / 3;      // K- 的数量
    double N_K_L = NPOT / 3;          // K_L 的数量

    // 计算每种介子的衰变分支比
    double Br_K_plus = Br_K_to_πφ("K+", θ, p_phi_CM_value);
    double Br_K_minus = Br_K_to_πφ("K-", θ, p_phi_CM_value);
    double Br_K_L = Br_K_to_πφ("K_L", θ, p_phi_CM_value);

    // 计算 N_phi
    return (NPOT / Nsim) * (N_K_plus * Br_K_plus + N_K_minus * Br_K_minus + N_K_L * Br_K_L) * ξ_K * Pphi;
}

// 从 CSV 文件中读取参数并计算 N_phi
void calculate_and_write_N_phi(const std::string& input_csv_file, const std::string& output_csv_file) {
    std::ifstream file(input_csv_file);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << input_csv_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::pair<double, double>> data;

    // 跳过标题行
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, ',')) {
            tokens.push_back(item);
        }

        if (tokens.size() == 2) {
            try {
                double m_phi = std::stod(tokens[0]);
                double θ_squared = std::stod(tokens[1]);
                data.emplace_back(m_phi, θ_squared);
            } catch (const std::exception& e) {
                std::cerr << "解析 CSV 文件中的数据时出错: " << e.what() << std::endl;
            }
        }
    }

    // 打开输出文件
    std::ofstream output_file(output_csv_file);
    if (!output_file.is_open()) {
        std::cerr << "无法打开输出文件: " << output_csv_file << std::endl;
        return;
    }

    // 写入标题
    output_file << "m_phi,θ_squared,N_phi\n";

    for (const auto& [m_phi, θ_squared] : data) {
        // 计算 N_phi
        double N_phi = calculate_N_phi(m_phi, θ_squared);

        // 写入到输出文件
        output_file << m_phi << "," << θ_squared << "," << N_phi << "\n";
    }

    std::cout << "计算完成，结果已写入到: " << output_csv_file << std::endl;
}

// 程序入口点
int main() {
    std::string input_csv_file = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/input/167PS191constraintsfilledarea.csv";
    std::string output_csv_file = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/167PS191output/167PS191constraintsfilledarea_Nphi.csv";

    calculate_and_write_N_phi(input_csv_file, output_csv_file);

    return 0;
}
