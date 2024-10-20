#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <iomanip>
#include <cmath>
#include <complex>
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/GM.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/Experiment_bound_167.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/Constants.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/FCNC.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/LTLauncher.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/SM.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/forward_scalar_rate.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/lorentz_vector.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/meson_decay.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/meson_rate.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/spdlog_wrapper.h"
#include "/home/fyq/DIYbyfyq/GM_FCNCDIY/include/utilities.h"

// 常量定义
const double NPOT = 2e6;

// 质量范围和实验 N_phi 区间结构体
struct MassRange {
    double m_phi_low;
    double m_phi_high;
    double N_phi_min;
    double N_phi_max;
};

// 从 CSV 文件中读取质量范围和实验 N_phi 区间数据
std::vector<MassRange> read_mass_ranges_from_csv(const std::string& file_path) {
    std::ifstream file(file_path);
    std::string line;
    std::vector<MassRange> ranges;

    // 跳过标题行
    std::getline(file, line);

    // 读取数据
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> tokens;

        while (std::getline(ss, item, ',')) {
            tokens.push_back(item);
        }

        if (tokens.size() == 4) {
            double m_phi_low = std::stod(tokens[0]);
            double m_phi_high = std::stod(tokens[1]);
            double N_phi_min = std::stod(tokens[2]);
            double N_phi_max = std::stod(tokens[3]);
            ranges.push_back({m_phi_low, m_phi_high, N_phi_min, N_phi_max});
        }
    }
    return ranges;
}

// 扫描质量和角度范围，并计算 N_phi
void scan_mass_and_theta(const std::vector<MassRange>& ranges, const std::string& output_file) {
    std::ofstream out_file(output_file);
    if (!out_file) {
        std::cerr << "无法打开文件 " << output_file << " 进行写入。" << std::endl;
        return;
    }

    out_file << "质量范围下限,质量范围上限,m_H5,theta_H,N_phi" << std::endl;
    out_file << std::fixed << std::setprecision(10);

    double N_K_plus = NPOT / 3;  
    double N_K_minus = NPOT / 3; 
    double N_K_L = NPOT / 3;     

    for (const auto& range : ranges) {
        double m_phi_low = range.m_phi_low;
        double m_phi_high = range.m_phi_high;
        double N_phi_min = range.N_phi_min;
        double N_phi_max = range.N_phi_max;

        for (double m_H5 = m_phi_low; m_H5 <= m_phi_high; m_H5 += (m_phi_high - m_phi_low) / 10) {
            for (double theta_H = -M_PI; theta_H <= M_PI; theta_H += 0.1) {
                // 使用 GM 计算 N_phi
                GM mod;
                mod.include_partial_two_loop(false);
                mod.set_input3(10, 10, m_H5, 1e-5, 0.01, 100, 0); // 这里假设的常量值

                int K_choice0 = 0;
                int K_choice1 = 1;
                double xi_H50_W = mod.get_H50_x_W();
                std::complex<double> xi_H50_ds = mod.get_H50_xij_down(FIRST, SECOND);
                std::complex<double> Br_Kpm = Br_K_PiH50(m_H5, xi_H50_W, xi_H50_ds, K_choice0);
                std::complex<double> Br_KL = Br_K_PiH50(m_H5, xi_H50_W, xi_H50_ds, K_choice1);
                std::complex<double> N_phi_GM = calculate_N_phi(Br_Kpm, Br_KL);

                if (N_phi_GM.real() >= N_phi_min && N_phi_GM.real() <= N_phi_max) {
                    out_file << m_phi_low << "," << m_phi_high << "," << m_H5 << "," << theta_H << "," << N_phi_GM.real() << std::endl;
                }
            }
        }
    }
}

int main() {
    // 输入文件和输出文件路径
    std::string input_file = "/home/fyq/Beam dump/CHARM/167/data/N_phi_output_167.csv";
    std::string output_file = "/home/fyq/Beam dump/CHARM/167/data/N_phi_scanned_results.csv";

    // 读取质量范围数据
    auto ranges = read_mass_ranges_from_csv(input_file);

    // 执行 GM 模型的计算和扫描
    scan_mass_and_theta(ranges, output_file);

    std::cout << "Data written to " << output_file << std::endl;
    return 0;
}
