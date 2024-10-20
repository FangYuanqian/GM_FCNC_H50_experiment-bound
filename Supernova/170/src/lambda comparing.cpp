#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector> // 加入vector头文件

// 常量定义
const double c = 3.0e10;          // cm/s, 光速
const double erg_to_MeV = 6.2415e5; // 1 erg = 6.2415e5 MeV
const double PI = 3.141592653589793;    // 圆周率
const double sqrt2 = sqrt(2.0);         // 根号2
const double m_pi = 135.0;        // MeV, pion mass
const double m_N = 940.0;         // MeV, nucleon mass
const double p_F = 200.0;         // MeV, Fermi momentum
const double T_SN = 30.0;         // MeV, supernova temperature
const double v_Higgs = 246.0 * 1000.0; // 转换为 MeV, Higgs vev

// 输入参数
double f_Tu_p = 0.023;
double f_Tu_n = 0.019;
double f_Td_p = 0.034;
double f_Td_n = 0.041;
double f_Ts_N = 0.14;

// 计算 G_phi(u)
double G_phi(double u) {
    double term1 = 1 - (5.0 / 2.0) * pow(u, 2) - (35.0 / 22.0) * pow(u, 4);
    double term2 = (5.0 / 64.0) * (28.0 * pow(u, 3) + 5.0 * pow(u, 5)) * atan(2.0 / u);
    double term3 = (5.0 / 64.0) * sqrt(2.0) * pow(u, 6) / sqrt(2.0 + pow(u, 2)) * atan(2.0 * sqrt(2.0 * (pow(u, 2) + 2)) / pow(u, 2));
    return term1 + term2 + term3;
    std::cout << "G_phi = " << term1 + term2 + term3 << std::endl;
}

// 计算 C_phiN
double C_phiN(double sin_theta, double f_Tu_p, double f_Tu_n, double f_Td_p, double f_Td_n, double f_Ts_N) {
    double f_Tu_N = f_Tu_p + f_Tu_n; // u夸克的贡献
    double f_Td_N = f_Td_p + f_Td_n; // d夸克的贡献
    double f_Tq_sum = f_Tu_N + f_Td_N + f_Ts_N; // 夸克总贡献
    double f_TG = 1.0 - f_Tq_sum; // 其他贡献

    double C_phiN_value = sqrt(2.0) * m_N / v_Higgs * ((6.0 / 27.0) * f_TG + f_Tq_sum) * sin_theta;

    return C_phiN_value;
}


// 计算 rho_phi
double calculate_rho_phi(double T_SN) {
    return (pow(PI, 2) * pow(T_SN, 4)) / 30.0;
}

// 计算 lambda_phi^{-1}
double calculate_lambda_phi_inverse(double Q_phi, double rho_phi) {
    return Q_phi / rho_phi * 5.06e10;
}

// 计算 lambda_phi
double calculate_lambda_phi1(double lambda_phi_inverse) {
    return 1.0 / lambda_phi_inverse;
}

// 计算 F
double calculate_F(double sin_theta) {
    return v_Higgs / sin_theta;
}

// 计算 lambda_phi2
double calculate_lambda_phi2(double F) {
    return 10.0 * 100.0 * pow(F / 1e9, 2);
}

int main() {
    std::vector<double> sin_theta_squared_values = {1}; // 替换为你需要计算的 sin(θ)^2 值

    for (double sin_theta_squared : sin_theta_squared_values) {
        std::cout << std::fixed << std::setprecision(10); // 设置输出精度
        std::cout << "Calculating for sin(θ)^2 = " << sin_theta_squared << std::endl;

        double sin_theta = sqrt(sin_theta_squared); // 计算 sin(theta)

        // 计算 G_phi(u)
        double u = m_pi / p_F;
        double G_phi_value = G_phi(u);

        // 计算 C_phiN
        double C_phiN_value = C_phiN(sin_theta, f_Tu_p, f_Tu_n, f_Td_p, f_Td_n, f_Ts_N);
        // 计算 Q_phi
        double Q_phi = pow(C_phiN_value, 2) * (11.0 / pow(15.0 * PI, 3)) * pow(T_SN / m_pi, 4) * pow(p_F, 5) * G_phi_value;

        // 输出调试信息
        std::cout << "C_phiN_value^2: " << pow(C_phiN_value, 2) << std::endl;
        std::cout << "Coefficient (330 / (15.0 * PI)^3): " << (330 / pow(15.0 * PI, 3) / pow(PI, 2)) << std::endl;
        std::cout << "(1 / m_pi)^4: " << pow(1 / m_pi, 4) << std::endl;
        std::cout << "(p_F)^5: " << pow(p_F, 5) << std::endl;
        std::cout << "111:" << pow(C_phiN_value, 2) * (330 / pow(15.0 * PI, 3) / pow(PI, 2)) * pow(1 / m_pi, 4) * pow(p_F, 5) * G_phi_value << std::endl;

        // 计算 Q_phi 的最终值
        std::cout << "Q_phi value: " << Q_phi << std::endl;

        // 计算 rho_phi
        double rho_phi = calculate_rho_phi(T_SN);
        std::cout << "rho_phi: " << rho_phi << std::endl;
        std::cout << "Q_phi / rho_phi : " << Q_phi / rho_phi  << std::endl;
        // 计算 lambda_phi^{-1}
        double lambda_phi_inverse1 = calculate_lambda_phi_inverse(Q_phi, rho_phi);
        std::cout << "lambda_phi^{-1}: " << Q_phi / rho_phi * 5.06e10 << std::endl;
        // 计算 λ_φ1
        double lambda_phi1 = calculate_lambda_phi1(lambda_phi_inverse1);
        std::cout << "lambda_phi: " << 1 / (Q_phi / rho_phi * 5.06e10) << std::endl;
        // 计算 F 和 λ_φ2
        double F = calculate_F(sin_theta);
        double lambda_phi2 = calculate_lambda_phi2(F);

        // 输出结果
        std::cout << "λ_φ1 (from Q_phi) for sin(θ)^2 = " << sin_theta_squared << ": " << lambda_phi1 << " cm" << std::endl;
        std::cout << "λ_φ2 (from F) for sin(θ)^2 = " << sin_theta_squared << ": " << lambda_phi2 << " cm" << std::endl;

        // 计算 λ_φ1 和 λ_φ2 的倍数关系
        if (lambda_phi1 != 0) { // 避免除以零
            double ratio = lambda_phi2 / lambda_phi1;
            std::cout << "Ratio (λ_φ2 / λ_φ1) for sin(θ)^2 = " << sin_theta_squared << ": " << ratio << std::endl;
        } else {
            std::cout << "λ_φ1 is zero, cannot calculate the ratio." << std::endl;
        }

        std::cout << "---------------------------------------" << std::endl;
    }

    return 0;
}
