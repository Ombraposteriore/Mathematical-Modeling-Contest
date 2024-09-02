import cv2 as cv
generation =151
import matplotlib.pyplot as plt

# 最终结果保存为视频
for iteration in range(1,generation,1):
    read_name = r'TSP_midResult\{0}.jpg'.format(iteration)
    img = cv.imread(read_name)
    plt.imshow(img)
    fps = 1
    videoWrite = cv.VideoWriter('MySaveVideo.avi', cv.VideoWriter_fourcc('I', '4', '2', '0'),
                                    fps,
                                    (640, 480), )

    videoWrite.write(img)
    print(read_name)