import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt




# class Encode():

# 生成编码，传入的参数为num_job, num_machine
def createChrome(num_job, num_machine):
    a = []
    for i in range(1, num_job + 1):
        a.append(i)

    chrome = []
    for j in range(num_machine):
        chrome += a  # 生成初始的编码
    random.shuffle(chrome)  # 将编码乱序排列，生成作业码的乱序
    return chrome

    # 根据种群数量和作业数量和设备数量生成初始种群


def createChromes(num_job, num_machine, pop):
    chromes = np.zeros((pop, num_job * num_machine), dtype=int)
    for i in range(pop):
        chromes[i, :] = createChrome(num_job, num_machine)
    return chromes


# 初始化工时设备数据
def initData(data):
   # m = data
    for i in range(0, data.shape[1], 2):
        data[:, i] = data[:, i] + 1
    return data


def createSchedule(data, chrome):
    num_job = np.size(data, 0)
    num_machine = np.size(data, 1) # 2
    schedule = np.zeros((num_job * num_machine, 5))

    # 定义中间数组
    jobCanStartTime = np.zeros((1, num_job), dtype=int)
    jobProcessId = np.ones((1, num_job), dtype=int)

    # 获取编码中当前的作业号
    chrome = createChrome(num_job, num_machine)
    data = initData(data)
    for i in range(num_job * num_machine):
        nowJobId = chrome[i] - 1
        nowProcessId = jobProcessId[0, nowJobId]
        nowMachId = data[nowJobId, 2 * nowProcessId - 2]
        nowProcTime = data[nowJobId, 2 * nowProcessId - 1]

        machSch = schedule[schedule[:, 1] == nowMachId, :]
        jobCanST = jobCanStartTime[0, nowJobId]
        if np.size(machSch, 0) == 0:  # 该工件还未安排作业
            startTime = jobCanStartTime[0, nowJobId]
            endTime = startTime + nowProcTime
        else:  # 设备已安排了工作
            machSch = machSch[np.argsort(machSch[:, 3]), :]
            rows = np.size(machSch, 0)
            # 处理第一行已排作业，检查是否能将当前作业排到之前
            done = 0
            if jobCanST < machSch[0, 3]:
                if machSch[0, 3] - jobCanST > nowProcTime:
                    startTime = jobCanST
                    endTime = startTime + nowProcTime
                    done = 1
            if done == 0:

                for j in range(rows):
                    if jobCanStartTime[0, nowJobId] < machSch[j, 3]:
                        if machSch[j, 3] - max(jobCanST, machSch[j - 1, 4]) > nowProcTime:
                            startTime = max(jobCanST, machSch[j - 1, 4])
                            endTime = startTime + nowProcTime
                            done = 1
                            break
            if done == 0:  # 表示该作业不能排到该设备已有作业之前
                startTime = max(jobCanST, machSch[rows - 1, 4])
                endTime = startTime + nowProcTime

        schedule[i, 0] = nowJobId + 1
        schedule[i, 1] = nowMachId
        schedule[i, 2] = nowProcessId
        schedule[i, 3] = startTime
        schedule[i, 4] = endTime
        jobCanStartTime[0, nowJobId] = endTime
        jobProcessId[0, nowJobId] += 1

    return schedule


# 根据调度方案绘制甘特图
def drawGant(schedule):
    rows = np.size(schedule, 0)
    num_job = int(max(schedule[:, 1]))
    mycolor = np.random.random((num_job, 3))
    mycolor = list(mycolor)

    for i in range(rows):
        x = [schedule[i, 3], schedule[i, 4]]
        y = [schedule[i, 1], schedule[i, 1]]
        n = int(schedule[i, 0])
        plt.plot(x, y, linewidth=8.0, color=mycolor[n - 1])
    plt.show()


num_job = 5 # 工件数目
num_machine = 4  #机器数目
pop = 5 #种群数

# data1 = [
#     [3, 1, 2, 3, 1, 6, 3, 7, 5, 3, 4, 6],
#     [1, 8, 2, 5, 4, 10, 5, 10, 0, 10, 3, 4],
#     [2, 5, 3, 4, 5, 8, 0, 9, 1, 1, 4, 7],
#     [1, 5, 0, 5, 2, 5, 3, 3, 4, 8, 5, 9],
#     [2, 9, 1, 3, 4, 5, 5, 4, 0, 3, 3, 1],
#     [1, 3, 3, 3, 5, 9, 0, 10, 4, 4, 2, 1],
# ]

data1 = [
    [31,41,25,30],
    [19,55,3,34],
    [23,42,27,6],
    [13,22,14,13],
    [33,5,57,19]
]
data = np.array(data1)
e = createChromes(num_job, num_machine, pop)
s = createSchedule(data, e)

print(s)#打印调度顺序
d = drawGant(s)#打印甘特图