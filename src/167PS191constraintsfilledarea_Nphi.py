import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取 C++ 生成的 CSV 文件
data_file = "/home/fyq/DIYbyfyq/H5/data/167PS191output/167PS191constraintsfilledarea_Nphi.csv"
df = pd.read_csv(data_file)

# 创建一个散点图
plt.figure(figsize=(10, 6))

# 计算 N_phi 的最小值和最大值
min_N_phi = df['N_phi'].min()
max_N_phi = df['N_phi'].max()

# 检查数据范围
print(f"Minimum N_phi: {min_N_phi}")
print(f"Maximum N_phi: {max_N_phi}")

# 绘制散点图，确保 N_phi 用于颜色映射
scatter = plt.scatter(df['m_phi'], df['θ_squared'], c=df['N_phi'], cmap='plasma', s=1, alpha=0.6, vmin=min_N_phi, vmax=max_N_phi)

# 添加颜色条
cbar = plt.colorbar(scatter, label="N_phi")
cbar.ax.tick_params(labelsize=8)  # 调整颜色条字体大小

# 设置颜色条的刻度为对数刻度，并确保覆盖到 N_phi 的范围
cbar.set_ticks(np.logspace(np.log10(min_N_phi), np.log10(max_N_phi), num=2))  # 设定对数刻度
cbar.set_ticklabels([f'{tick:.1e}' for tick in np.logspace(np.log10(min_N_phi), np.log10(max_N_phi), num=2)])  # 格式化刻度标签为科学计数法

# 设置坐标轴为对数坐标
plt.xscale('log')
plt.yscale('log')

# 设置坐标轴刻度，手动指定避免重叠
plt.xticks([3e-2, 0.1, 1], ['0.003', '0.1', '1']) 
plt.yticks([1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3], 
           ['1e-10', '1e-9', '1e-8', '1e-7', '1e-6', '1e-5', '1e-4', '1e-3'], rotation=45)  # 旋转刻度标签以减少重叠

# 设置图表标题和标签
plt.title('167PS191_N_phi Distribution Plot')
plt.xlabel('m_phi (GeV)')
plt.ylabel('θ_squared')

# 保存图像或显示图表
plt.savefig("/home/fyq/DIYbyfyq/H5/data/167PS191output/167PS191constraintsfilledarea_Nphi_scatter.png")  # 保存为图片
plt.show()  # 显示图表
