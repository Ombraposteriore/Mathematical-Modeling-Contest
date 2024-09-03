# import pandas as pd
# import matplotlib.pyplot as plt
#
# # 设置中文字体
# plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
# plt.rcParams['axes.unicode_minus'] = False  # 解决负号无法显示的问题
#
# # 读取Excel文件（假设文件名为data.xlsx，并且数据在第一个工作表）
# df = pd.read_excel('数据集.xlsx')  # 修改为你的文件名
#
# # 提取第二列数据（假设是数据列）
# data = df.iloc[:, 1]  # 第二列是数据
#
# # 绘制图表，调整图的大小
# plt.figure(figsize=(30, 20))  # (x,y) 宽x，高y
# plt.plot(data, marker='o', linestyle='-')
# plt.title('数据按顺序绘制')
# plt.xlabel('索引')
# plt.ylabel('数据值')
# plt.grid(True)
# plt.show()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from keras import Sequential
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from tensorflow.keras.layers import Dense
from tcn import TCN  # pip install keras-tcn

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体
plt.rcParams['axes.unicode_minus'] = False  # 解决负号无法显示的问题

# 读取Excel文件（假设文件名为数据集.xlsx，并且数据在第一个工作表）
df = pd.read_excel('数据集.xlsx')  # 修改为你的文件名

# 提取第二列数据（假设是数据列）
data = df.iloc[:, 1].values  # 第二列是数据

# 数据预处理：创建训练和测试数据集
def create_sequences(data, seq_length):
    X = []
    y = []
    for i in range(len(data) - seq_length):
        X.append(data[i:i + seq_length])
        y.append(data[i + seq_length])
    return np.array(X), np.array(y)

# 设置序列长度
seq_length = 10  # 可以调整序列长度

# 创建时间序列样本数据
X, y = create_sequences(data, seq_length)

# 将数据形状调整为模型输入的格式 (样本数, 时间步长, 特征数)
X = X.reshape((X.shape[0], X.shape[1], 1))

# 划分训练集、验证集和测试集（80% 训练，10% 验证，10% 测试）
train_size = int(len(X) * 0.8)
val_size = int(len(X) * 0.1)

X_train, y_train = X[:train_size], y[:train_size]
X_val, y_val = X[train_size:train_size + val_size], y[train_size:train_size + val_size]
X_test, y_test = X[train_size + val_size:], y[train_size + val_size:]

# 构建 TCN 模型
def create_tcn_model(input_shape):
    model = Sequential([
        TCN(input_shape=input_shape, return_sequences=False),  # TCN层
        Dense(1, activation='linear')  # 输出层
    ])
    model.compile(optimizer='adam', loss='mse')
    return model

# 创建并训练 TCN 模型
model = create_tcn_model(input_shape=(seq_length, 1))
model.summary()
model.fit(X_train, y_train, epochs=20, batch_size=32, validation_data=(X_val, y_val))

# 使用测试集进行预测
y_pred = model.predict(X_test)

# 计算评估指标
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_test, y_pred)
def mean_absolute_percentage_error(y_true, y_pred):
    return np.mean(np.abs((y_true - y_pred) / y_true)) * 100
mape = mean_absolute_percentage_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f'MSE: {mse:.4f}')
print(f'RMSE: {rmse:.4f}')
print(f'MAE: {mae:.4f}')
print(f'MAPE: {mape:.2f}%')
print(f'R^2: {r2:.4f}')

# 可视化预测结果
plt.figure(figsize=(12, 6))

# 绘制原始数据的折线图
plt.plot(range(len(data)), data, label='原始数据', color='blue')

# 绘制预测结果的折线图
plt.plot(range(len(data) - len(y_pred), len(data)), y_pred.flatten(), label='预测值', color='red', linestyle='-')

plt.title('时间序列数据预测（TCN）')
plt.xlabel('时间步')
plt.ylabel('数据值')
plt.legend()
plt.grid(True)
plt.show()

'''
MSE 和 RMSE 对大的误差更敏感，适用于需要惩罚大误差的场景。
MAE 不容易受异常值的影响，适用于对异常值不敏感的场景。
MAPE 是一个无量纲的指标，适用于不同尺度数据的比较，但在实际值接近零时会出现问题。
R^2值 表示拟合优度，适用于判断模型的解释能力。

结果分析：
结果分析
1.MSE (均方误差): 0.0008
MSE 是衡量模型预测值与真实值之间的平均平方误差。越接近 0，模型预测越准确。0.0008 是一个相对较小的值，这表明模型的预测值与真实值之间的误差较小。

2.RMSE (均方根误差): 0.0282
RMSE 是 MSE 的平方根，表示预测误差的标准差。RMSE 越接近 0，模型越好。0.0282 也是一个相对较低的值，说明模型预测误差相对较小，预测的值和真实值之间非常接近。

3.MAE (平均绝对误差): 0.0234
MAE 是预测值与真实值之间绝对差的平均值。它是一个比 MSE 和 RMSE 更具解释性的指标，表示每个预测的平均绝对误差。在 0.0234 这个范围内，模型误差也是相对较小的。

4.MAPE (平均绝对百分比误差): 11.30%
MAPE 反映了预测值与真实值的相对误差。11.30% 表示预测值与真实值之间的误差平均为 11.3%。这个值对于许多实际应用来说是可以接受的，但在某些对预测精度要求更高的场景中，可能还需要进一步优化。

5.R^2(决定系数): 0.7657
R^2 是模型拟合优度的衡量标准，取值范围为 0 到 1。0.7657 表示模型能够解释约 76.57% 的数据方差，这说明模型有一定的预测能力，但并不完美。
通常，R^2值在 0.7 到 0.9 之间表示模型有良好的拟合能力。对于时间序列数据，这样的R^2值是相对合理的，但也可能存在进一步优化的空间。
'''