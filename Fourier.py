import numpy as np
import cv2
import sys
from matplotlib import pyplot as plt
import csv

#def nothing(x):
#    pass

#print(cv2.__version__) # print the version Opencv (=3.4.2)
file_src = '/Users/konrai/Dropbox/Konrai_exp/konara_mizunara/Threshold/449-konara-img023_1.jpg'
#file_dst ='/Users/konrai/Dropbox/My Mac (Konrai.local)/Desktop/dst.jpg'

#img_src = cv2.imread(file_src, 1) # load as color image
img_src = cv2.imread(file_src, 0) # load as gray scale image
# 読み込んだ画像の高さと幅を取得
height = img_src.shape[0]
width = img_src.shape[1]
#print('height:', height)
#print('width:', width)

#img_src = cv2.resize(img_src, (round(width/2), round(height/2))) # Resize

cv2.namedWindow('src', cv2.WINDOW_NORMAL)
#cv2.namedWindow('dst', cv2.WINDOW_NORMAL)

# show input image
#cv2.imshow('src', img_src) # show imput image
#cv2.waitKey(0) # Wait until any key input
#cv2.destroyAllWindows()
#cv2.createTrackbar('threshold', 'src', 0, 255, nothing) # trackbar to threshold

#plt.imshow(img_src, cmap='gray', vmin=0, vmax=255)
#plt.colorbar()
#plt.show()

### Closing (Removing noises) うまくいかん。。。
#element8 = np.array([[1,1,1], [1,1,1], [1,1,1]], np.uint8) # 8-neighborhood (8近傍)
#img_tmp = cv2.morphologyEx(img_src, cv2.MORPH_CLOSE, element8)
#plt.imshow(img_tmp, cmap='gray', vmin=0, vmax=255)
#plt.colorbar()

# Rabelling
nlabel, img_lab = cv2.connectedComponents(image = img_src)
print('Number of labels:', nlabel)
for n in range(1, 10):
    print(n)
    img_tmp = cv2.compare(img_lab, n, cv2.CMP_EQ) # Extract image which is label n
    contours = cv2.findContours(img_tmp, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[1]
    area = cv2.contourArea(contours[0]) # area
    perimeter = cv2.arcLength(np.array(contours[0]), True)
    if area != 0 and perimeter != 0:
        roundness = 4 * np.pi * area / perimeter / perimeter
    else:
        roundness = 0
    print('Area:', area)
    print('Perimeter:', perimeter)
    print('Roundness:', roundness)
    plt.imshow(img_tmp, cmap='gray', vmin=0, vmax=255)
    plt.show()
#img_tmp = cv2.compare(img_lab, 1, cv2.CMP_EQ)
#plt.imshow(img_lab, cmap='gray', vmin=0, vmax=255)
#plt.colorbar()
#plt.show()

#while(1):
#    cv2.imshow('src', img_src)
#    k = cv2.waitKey(1) & 0xFF
#    if k == 27:
#        break

    # get current positions of four trackbars
#    thresh = cv2.getTrackbarPos('threshold', 'src')
#    ret, img_src = cv2.threshold(img_src, thresh, 255, cv2.THRESH_BINARY)

#cv2.destroyAllWindows()

# Binaryzation 2値化
#thresh = 100 # threshold
#ret, img_dst = cv2.threshold(img_src, thresh, 255, cv2.THRESH_BINARY)

#cv2.imshow('dst', img_dst) # show output image
#cv2.waitKey(0) # Wait until any key input
#cv2.destroyAllWindows()

#plt.imshow(img_dst, cmap='gray', vmin=0, vmax=255)
#plt.colorbar()
#plt.show()

### Generate example image (circle)
#img_src = np.zeros((300, 500, 1), dtype=np.uint8)
#height = img_src.shape[0] # image height
#width = img_src.shape[1] # image width
#print('height: ', height)
#print('width: ', width)
#cv2.circle(img_src, (round(width/2), round(height/2)), 50, color=255, thickness=-1) # add circle on img_src

### Get some properties
#x, y, w, h = cv2.boundingRect(img_src) # 外接長方形を求める。Rect型で(x, y)：長方形の左上頂点座標, width, heightはそれぞれ長方形の横幅、縦幅
#aspectratio = float(h) / w # height / width
#print(aspectratio)
#cv2.rectangle(img_src, (x,y), (x+w, y+h), 128) # show bounding rectangle

### Find contours (輪郭の取得)
image, contours, hierarchy =  cv2.findContours(image=img_src, mode=cv2.RETR_EXTERNAL, method=cv2.CHAIN_APPROX_NONE) # The most external contour is only stored
#print('image: ', image, type(image))
#print('contours: ', contours, type(contours))
#print('hierarchy: ', hierarchy, type(hierarchy))

#print(contours[0][1][0])
#print(len(contours[0][:]))
l =list() # generate list in which contour coorinate(x, y) is added
#print(l)
for i in range(len(contours[0][:])):
    print(contours[0][i][0])
    l.append(contours[0][i][0])
with open('/Users/konrai/Dropbox/My Mac (Konrai.local)/Desktop/contours.csv', 'w') as f:
    writer = csv.writer(f) # save as csv file
    writer.writerows(l)

#np.savetxt('/Users/konrai/Dropbox/My Mac (Konrai.local)/Desktop/contours.csv', contours[0][0:5], delimiter=',')

#plt.imshow(img_src, cmap='gray', vmin=0, vmax=255) # show the input (original) image
#plt.colorbar() # add color bar
#plt.show()

# draw contours (輪郭の描画)
img_dst = cv2.drawContours(image=img_src, contours=contours, contourIdx=-1, color=128) # show out put(contour) image as gray line
plt.imshow(img_dst, cmap='gray', vmin=0, vmax=255)
plt.colorbar()
#plt.show()
