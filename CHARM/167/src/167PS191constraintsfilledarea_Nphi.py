import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取 CSV 文件
data_file = "/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/167PS191output/167PS191constraintsfilledarea_Nphi.csv"
df = pd.read_csv(data_file)

# 创建一个散点图
plt.figure(figsize=(10, 6))

# 计算 N_phi 的最小值和最大值
min_N_phi = df['N_phi'].min()
max_N_phi = df['N_phi'].max()

# 绘制散点图，确保 N_phi 用于颜色映射
scatter = plt.scatter(df['m_phi'], df['θ_squared'], c=df['N_phi'], cmap='plasma', s=1, alpha=0.6, vmin=min_N_phi, vmax=max_N_phi)

# 添加颜色条
cbar = plt.colorbar(scatter, label="N_phi")
cbar.ax.tick_params(labelsize=8)

# 设置颜色条的刻度为对数刻度，并显示更多的刻度
log_ticks = np.logspace(np.log10(min_N_phi), np.log10(max_N_phi), num=3)  # 增加更多刻度
cbar.set_ticks(log_ticks)
cbar.set_ticklabels([f'{tick:.1e}' for tick in log_ticks])

# 设置坐标轴为对数坐标
plt.xscale('log')
plt.yscale('log')

# 设置坐标轴刻度
plt.xticks([3e-2, 0.1, 1], ['0.003', '0.1', '1'])
plt.yticks([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3], 
           ['1e-10', '1e-9', '1e-8', '1e-7', '1e-6', '1e-5', '1e-4', '1e-3'], rotation=45)

# 设置图表标题和标签
plt.title('167PS191_N_phi Distribution Plot')
plt.xlabel('m_phi (GeV)')
plt.ylabel('θ_squared')

# 为每个数量级选择一个点标注 N_phi 数值
log_N_phi_levels = np.arange(np.floor(np.log10(min_N_phi)), np.ceil(np.log10(max_N_phi)) + 1)

for log_level in log_N_phi_levels:
    target_value = 10 ** log_level
    closest_point = df.iloc[(df['N_phi'] - target_value).abs().argmin()]
    
    m_phi = closest_point['m_phi']
    theta_sq = closest_point['θ_squared']
    N_phi = closest_point['N_phi']
    
    plt.text(m_phi, theta_sq, f'{N_phi:.1e}', fontsize=8, color='red', ha='center', va='bottom')

# 保存图像或显示图表
plt.savefig("/home/fyq/DIYbyfyq/Beam dump/CHARM/167/data/167PS191output/167PS191constraintsfilledarea_Nphi_scatter")
plt.show()
