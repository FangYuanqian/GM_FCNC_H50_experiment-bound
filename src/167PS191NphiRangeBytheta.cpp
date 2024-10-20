#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>
#include <utility>

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
        m_K = 0.497614;  // K_L 的质量，单位：GeV/c^2
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

// 从 CSV 文件中读取数据
std::vector<std::pair<double, double>> read_m_phi_theta_squared_from_csv(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::vector<std::pair<double, double>> data;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, ',')) {
            tokens.push_back(item);
        }

        if (tokens.size() == 2) {
            double m_phi = std::stod(tokens[0]);
            double theta_squared = std::stod(tokens[1]);
            data.emplace_back(m_phi, theta_squared);
        }
    }
    return data;
}

// 计算每个 θ^2 区间内的 N_phi 的上下边缘，并输出到终端和 CSV
void calculate_N_phi_in_ranges(const std::vector<std::string>& csv_files, const std::string& output_file) {
    std::ofstream out_file(output_file);
    if (!out_file) {
        std::cerr << "Unable to open file " << output_file << " for writing." << std::endl;
        return;
    }

    // 输出格式：θ^2 Lower Bound, θ^2 Upper Bound, N_phi Min, N_phi Max
    out_file << "θ^2 Lower Bound,θ^2 Upper Bound,N_phi Min,N_phi Max" << std::endl;

    for (const auto& csv_file : csv_files) {
        auto data = read_m_phi_theta_squared_from_csv(csv_file);

        double min_theta_squared = std::numeric_limits<double>::max();
        double max_theta_squared = std::numeric_limits<double>::lowest();

        for (const auto& [_, θ_squared] : data) {
            min_theta_squared = std::min(min_theta_squared, θ_squared);
            max_theta_squared = std::max(max_theta_squared, θ_squared);
        }

        int num_points = 10; 
        double step = (max_theta_squared - min_theta_squared) / (num_points - 1);

        std::cout << "File: " << csv_file << " - N_phi range:" << std::endl;

        for (int i = 0; i < num_points - 1; ++i) {
            double theta_squared_low = min_theta_squared + i * step;
            double theta_squared_high = min_theta_squared + (i + 1) * step;

            double min_N_phi = std::numeric_limits<double>::max();
            double max_N_phi = std::numeric_limits<double>::lowest();

            for (const auto& [m_phi, θ_squared] : data) {
                if (θ_squared >= theta_squared_low && θ_squared <= theta_squared_high) {
                    double N_phi = calculate_N_phi(m_phi, θ_squared);

                    if (!std::isnan(N_phi) && !std::isinf(N_phi)) {
                        min_N_phi = std::min(min_N_phi, N_phi);
                        max_N_phi = std::max(max_N_phi, N_phi);
                    }
                }
            }

            // 终端输出
            std::cout << "[" << theta_squared_low << ", " << theta_squared_high
                      << "] -> N_phi Range: [" << min_N_phi << ", " << max_N_phi << "]" << std::endl;

            // 输出到 CSV 文件
            out_file << theta_squared_low << "," << theta_squared_high
                     << "," << min_N_phi << "," << max_N_phi << std::endl;
        }
    }

    out_file.close();
}

// 程序入口点
int main() {
    std::vector<std::string> csv_files = {
        "/home/fyq/DIYbyfyq/H5/data/input/167PS191constraintsupperline.csv",
        "/home/fyq/DIYbyfyq/H5/data/input/167PS191constraintsdownline.csv"
    };
    std::string output_file = "/home/fyq/DIYbyfyq/H5/data/167PS191output/167PS191NphiRangeBytheta_output.csv";
    calculate_N_phi_in_ranges(csv_files, output_file);
    return 0;
}
